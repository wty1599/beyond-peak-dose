\set ON_ERROR_STOP on

BEGIN;

CREATE SCHEMA IF NOT EXISTS project_d4;

-- Original study definition: latest hospital discharge across all admissions.
-- This is not the latest observation across all clinical source tables.
DROP TABLE IF EXISTS project_d4.d4_subject_last_record;
CREATE TABLE project_d4.d4_subject_last_record AS
SELECT
    a.subject_id,
    MAX(a.dischtime) AS last_documented_dischtime
FROM mimiciv_hosp.admissions AS a
GROUP BY a.subject_id;
CREATE UNIQUE INDEX d4_subject_last_record_idx
    ON project_d4.d4_subject_last_record (subject_id);

COMMIT;
