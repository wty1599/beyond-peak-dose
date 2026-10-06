\set ON_ERROR_STOP on
\pset pager off
\pset null '[NULL]'

BEGIN;

CREATE SCHEMA IF NOT EXISTS project_d2v2;

-- 1. Blood lactate only. LDH and non-blood fluids are retained as excluded
-- candidates in the export so the label-based resolution is auditable.
DROP TABLE IF EXISTS project_d2v2.d2v2_map_lactate;
CREATE TABLE project_d2v2.d2v2_map_lactate AS
SELECT
    itemid,
    label,
    fluid,
    category,
    CASE
        WHEN lower(label) = 'lactate'
         AND lower(fluid) = 'blood'
        THEN true ELSE false
    END AS include_main,
    CASE
        WHEN lower(label) = 'lactate'
         AND lower(fluid) = 'blood'
        THEN 'BLOOD_LACTATE'
        ELSE 'EXCLUDE_NON_BLOOD_OR_LDH'
    END AS mapping_reason
FROM mimiciv_hosp.d_labitems
WHERE lower(label) ~ 'lactate|lactic acid'
ORDER BY itemid;

-- 2. Positive inotropes used for the cardiogenic-shock exclusion.
DROP TABLE IF EXISTS project_d2v2.d2v2_map_inotrope;
CREATE TABLE project_d2v2.d2v2_map_inotrope AS
SELECT
    itemid,
    label,
    abbreviation,
    category,
    linksto,
    CASE
        WHEN lower(coalesce(label, '') || ' ' || coalesce(abbreviation, ''))
             ~ 'dobutamine|milrinone'
         AND linksto = 'inputevents'
        THEN true ELSE false
    END AS include_main,
    CASE
        WHEN lower(coalesce(label, '') || ' ' || coalesce(abbreviation, ''))
             ~ 'dobutamine' THEN 'DOBUTAMINE'
        WHEN lower(coalesce(label, '') || ' ' || coalesce(abbreviation, ''))
             ~ 'milrinone' THEN 'MILRINONE'
        ELSE 'EXCLUDE'
    END AS drug_class
FROM mimiciv_icu.d_items
WHERE lower(coalesce(label, '') || ' ' || coalesce(abbreviation, ''))
      ~ 'dobutamine|milrinone'
ORDER BY itemid;

-- 3. Any temporary or durable mechanical circulatory support device.
-- Generic "Assist Device" is deliberately excluded because it records
-- mobility aids/braces rather than MCS. Both chart and procedure sources are
-- retained, including device-specific activity parameters.
DROP TABLE IF EXISTS project_d2v2.d2v2_map_mcs;
CREATE TABLE project_d2v2.d2v2_map_mcs AS
WITH candidates AS (
    SELECT
        itemid,
        linksto,
        label,
        abbreviation,
        category,
        unitname,
        lower(
            coalesce(label, '') || ' ' ||
            coalesce(abbreviation, '') || ' ' ||
            coalesce(category, '')
        ) AS txt
    FROM mimiciv_icu.d_items
    WHERE linksto IN ('chartevents', 'procedureevents')
      AND lower(
            coalesce(label, '') || ' ' ||
            coalesce(abbreviation, '') || ' ' ||
            coalesce(category, '')
          ) ~
          'ecmo|oxygenator|sweep|impella|centrimag|tandem.?heart|intra.?aortic|iabp|heartmate|hm ii|ventricular assist|(^|[^a-z])vad([^a-z]|$)|durable vad|assist device'
)
SELECT
    itemid,
    linksto,
    label,
    abbreviation,
    category,
    unitname,
    CASE
        WHEN txt ~ 'tandem.?heart' THEN 'TANDEMHEART'
        WHEN txt ~ 'ecmo|oxygenator|sweep' THEN 'ECMO'
        WHEN txt ~ 'impella' THEN 'IMPELLA'
        WHEN txt ~ 'intra.?aortic|iabp' THEN 'IABP'
        WHEN txt ~ 'centrimag' THEN 'CENTRIMAG'
        WHEN txt ~ 'heartmate|hm ii|ventricular assist|(^|[^a-z])vad([^a-z]|$)|durable vad'
            THEN 'VAD_HEARTMATE'
        ELSE 'GENERIC_ASSIST_DEVICE'
    END AS mcs_type,
    CASE
        WHEN lower(label) = 'assist device' THEN false
        WHEN txt ~
             'ecmo|oxygenator|sweep|impella|centrimag|tandem.?heart|intra.?aortic|iabp|heartmate|hm ii|ventricular assist|(^|[^a-z])vad([^a-z]|$)|durable vad'
        THEN true ELSE false
    END AS include_main,
    CASE
        WHEN lower(label) = 'assist device'
            THEN 'EXCLUDE_GENERIC_MOBILITY_AID'
        WHEN linksto = 'procedureevents'
            THEN 'DEVICE_PROCEDURE_INTERVAL'
        ELSE 'DEVICE_SPECIFIC_CHART_ACTIVITY'
    END AS mapping_reason
FROM candidates
ORDER BY linksto, mcs_type, itemid;

-- 4. Cardiac critical-care units.
DROP TABLE IF EXISTS project_d2v2.d2v2_map_cardiac_unit;
CREATE TABLE project_d2v2.d2v2_map_cardiac_unit AS
SELECT
    first_careunit,
    CASE
        WHEN lower(first_careunit) ~
             'cardiac vascular intensive care|coronary care'
        THEN true ELSE false
    END AS include_cardiogenic_exclusion,
    CASE
        WHEN lower(first_careunit) ~ 'cardiac vascular intensive care'
            THEN 'CARDIAC_SURGICAL_ICU'
        WHEN lower(first_careunit) ~ 'coronary care'
            THEN 'CARDIAC_MEDICAL_ICU'
        ELSE 'NON_CARDIAC_FIRST_ICU'
    END AS unit_class
FROM (
    SELECT DISTINCT first_careunit
    FROM mimiciv_icu.icustays
) AS u
ORDER BY first_careunit;

-- 5. Shock/hemorrhage ICD candidates. Shock-specific codes are included.
-- Acute post-haemorrhagic anaemia or unspecified haemorrhage alone is not
-- considered sufficient to label the vasopressor episode hemorrhagic.
DROP TABLE IF EXISTS project_d2v2.d2v2_map_shock_icd;
CREATE TABLE project_d2v2.d2v2_map_shock_icd AS
WITH candidates AS (
    SELECT
        d.icd_code,
        d.icd_version,
        d.long_title,
        count(di.*)::bigint AS diagnosis_records
    FROM mimiciv_hosp.d_icd_diagnoses AS d
    LEFT JOIN mimiciv_hosp.diagnoses_icd AS di
      ON d.icd_code = di.icd_code
     AND d.icd_version = di.icd_version
    WHERE lower(d.long_title) ~
          'cardiogenic shock|hypovolemic shock|traumatic shock|hemorrhagic shock|acute posthemorrhagic anemia|haemorrhage, unspecified|hemorrhage, unspecified|hemorrhage, not elsewhere classified'
    GROUP BY d.icd_code, d.icd_version, d.long_title
)
SELECT
    icd_code,
    icd_version,
    long_title,
    diagnosis_records,
    CASE
        WHEN lower(long_title) ~ 'cardiogenic shock'
            THEN 'CARDIOGENIC_SHOCK'
        WHEN lower(long_title) ~
             'hypovolemic shock|traumatic shock|hemorrhagic shock'
            THEN 'HEMORRHAGIC_OR_HYPOVOLEMIC_SHOCK'
        ELSE 'HEMORRHAGE_NONSPECIFIC'
    END AS shock_class,
    CASE
        WHEN lower(long_title) ~
             'cardiogenic shock|hypovolemic shock|traumatic shock|hemorrhagic shock'
        THEN true ELSE false
    END AS include_main,
    CASE
        WHEN lower(long_title) ~
             'cardiogenic shock|hypovolemic shock|traumatic shock|hemorrhagic shock'
            THEN 'SHOCK_SPECIFIC_CODE'
        ELSE 'EXCLUDE_NOT_SHOCK_SPECIFIC'
    END AS mapping_reason
FROM candidates
ORDER BY shock_class, diagnosis_records DESC, icd_version, icd_code;

-- 6. The six approved crystalloid/colloid items, re-resolved from labels.
-- Dextran 70 remains visible as an inactive/zero-record exclusion rather than
-- being silently deleted.
DROP TABLE IF EXISTS project_d2v2.d2v2_map_fluid;
CREATE TABLE project_d2v2.d2v2_map_fluid AS
WITH candidates AS (
    SELECT
        d.itemid,
        d.label,
        d.abbreviation,
        d.category,
        d.unitname,
        d.linksto,
        count(i.*)::bigint AS raw_records,
        count(i.*) FILTER (
            WHERE i.amount > 0 AND lower(coalesce(i.amountuom, '')) = 'ml'
        )::bigint AS positive_ml_records
    FROM mimiciv_icu.d_items AS d
    LEFT JOIN mimiciv_icu.inputevents AS i
      ON d.itemid = i.itemid
    WHERE d.linksto = 'inputevents'
      AND lower(coalesce(d.label, '') || ' ' || coalesce(d.abbreviation, ''))
          ~ '(^|[^a-z])lr([^a-z]|$)|nacl 0[,.]9%|albumin 5%|albumin 25%|hetastarch.*6%|dextran 40|dextran 70|plasmalyte|normosol|isolyte'
    GROUP BY
        d.itemid, d.label, d.abbreviation, d.category, d.unitname, d.linksto
)
SELECT
    *,
    CASE
        WHEN lower(label) IN ('lr', 'nacl 0.9%', 'albumin 5%', 'albumin 25%',
                              'hetastarch (hespan) 6%', 'dextran 40')
        THEN true ELSE false
    END AS include_main,
    CASE
        WHEN lower(label) IN ('lr', 'nacl 0.9%') THEN 'CRYSTALLOID'
        WHEN lower(label) IN ('albumin 5%', 'albumin 25%',
                              'hetastarch (hespan) 6%', 'dextran 40')
            THEN 'COLLOID'
        WHEN lower(label) = 'dextran 70' THEN 'ZERO_RECORD_RETAINED'
        WHEN lower(label) ~ 'plasmalyte|normosol|isolyte'
            THEN 'BALANCED_CRYSTALLOID_NO_ACTIVE_RECORD'
        ELSE 'EXCLUDE_NON_APPROVED_OR_INACTIVE'
    END AS fluid_class
FROM candidates
ORDER BY include_main DESC, fluid_class, itemid;

-- 7. PRBC. Source volumes are in mL; the screening script converts cumulative
-- volume to unit-equivalents using the observed single-unit mode/median of
-- approximately 350 mL. This operationalisation is reported as a limitation.
DROP TABLE IF EXISTS project_d2v2.d2v2_map_prbc;
CREATE TABLE project_d2v2.d2v2_map_prbc AS
SELECT
    d.itemid,
    d.label,
    d.abbreviation,
    d.category,
    d.unitname,
    d.linksto,
    count(i.*)::bigint AS raw_records,
    count(i.*) FILTER (
        WHERE i.amount > 0
          AND lower(coalesce(i.amountuom, '')) = 'ml'
    )::bigint AS positive_ml_records,
    percentile_cont(0.5) WITHIN GROUP (
        ORDER BY i.amount
    ) FILTER (
        WHERE i.amount > 0
          AND lower(coalesce(i.amountuom, '')) = 'ml'
    ) AS median_positive_ml,
    CASE
        WHEN lower(d.label) = 'packed red blood cells'
         AND d.linksto = 'inputevents'
        THEN true ELSE false
    END AS include_main,
    CASE
        WHEN lower(d.label) = 'packed red blood cells'
            THEN 'PRBC_VOLUME_ML'
        ELSE 'EXCLUDE_INACTIVE_LEGACY_ITEM'
    END AS mapping_reason
FROM mimiciv_icu.d_items AS d
LEFT JOIN mimiciv_icu.inputevents AS i
  ON d.itemid = i.itemid
WHERE d.linksto = 'inputevents'
  AND lower(d.label) IN ('packed red blood cells', 'packed red cells')
GROUP BY
    d.itemid, d.label, d.abbreviation, d.category, d.unitname, d.linksto
ORDER BY include_main DESC, d.itemid;


COMMIT;
