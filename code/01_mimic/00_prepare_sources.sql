\set ON_ERROR_STOP on

-- Prepare study-derived tables before the read-only candidate export.
-- Relative includes are resolved from this file, not the caller's directory.
\ir 00_source_mappings.sql
\ir 00_source_exclusions.sql
\ir 00_code_status.sql
\ir 00_subject_last_record.sql
