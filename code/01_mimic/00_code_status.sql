\set ON_ERROR_STOP on

CREATE SCHEMA IF NOT EXISTS project_step4;

DROP TABLE IF EXISTS project_step4.step4_poe_code_status_orders;
CREATE TABLE project_step4.step4_poe_code_status_orders AS
SELECT
    p.subject_id,
    p.hadm_id,
    p.poe_id,
    p.ordertime,
    CASE
        WHEN BOOL_OR(
            lower(d.field_value)
                ~ 'do not resuscitate|dnr/dni|dnar'
        ) THEN 'LIMITATION_PRESENT'
        WHEN BOOL_OR(
            lower(d.field_value) ~ 'full code|resuscitate'
        ) THEN 'FULL_CODE_DOCUMENTED'
        ELSE 'UNKNOWN'
    END AS code_status
FROM mimiciv_hosp.poe AS p
INNER JOIN mimiciv_hosp.poe_detail AS d
    ON p.poe_id = d.poe_id
WHERE p.order_type = 'General Care'
  AND p.order_subtype = 'Code status'
GROUP BY p.subject_id, p.hadm_id, p.poe_id, p.ordertime;

CREATE INDEX step4_poe_code_status_lookup_idx
    ON project_step4.step4_poe_code_status_orders
       (hadm_id, ordertime, poe_id);

