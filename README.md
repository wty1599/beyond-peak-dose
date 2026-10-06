# Beyond peak dose

Analysis code and disclosure-protected aggregate results for *Beyond peak dose:
the 72-hour transition after vasopressor withdrawal in septic shock* by Wang
and colleagues, prepared for Annals of Intensive Care (publication pending).

**Code and aggregate-results package.** The original MIMIC last-record table
definition is included, and the fixed-cohort predictor extractors no longer
require historical patient-level feature caches. The portable routes have not
been tested end to end with source data. No archival DOI has been assigned.

## Data access

Each database requires its own access application and applicable data-use
agreement. This repository contains no patient-level data or fitted/imputation
objects. Access to code does not confer access to the source databases.

- [MIMIC-IV v3.1](https://physionet.org/content/mimiciv/3.1/)
- [SICdb v1.0.8](https://physionet.org/content/sicdb/1.0.8/)
- [AmsterdamUMCdb](https://amsterdammedicaldatascience.nl),
  [v1.0.2 deposited dataset](https://doi.org/10.17026/DANS-22U-F8VD)

`results/aggregate/` contains publication-facing summaries. Small counts are
shown as `<10`; linked counts, denominators and percentages are withheld where
needed to prevent reconstruction from this package. Some tables therefore
retain only effect estimates or continuous summaries. Exact event-resolution
curve data are not distributed; fixed-time probabilities are rounded to 0.1
percentage point. Other retained estimates keep their original precision.

The PDFs are public-disclosure derivatives, not substitutes for the manuscript
figures: linked flow counts are withheld and time-course plots show only
released time points without interpolation. These restrictions are specific to
this repository; the manuscript remains the record of the complete study.

## Setup and execution

Use the repository root as the working directory. See
`environment/README.md` and `environment/sessionInfo.txt` for dependencies and
recorded versions. Copy `config/paths.example.yml` to `config/paths.yml`, then
replace the source-location placeholders. The file uses JSON-compatible YAML
so that both languages read the same values without a second YAML dependency.
Set `BEYOND_PATHS_CONFIG` only if using another local configuration file.

Set PostgreSQL credentials locally through `MIMIC_PG_HOST`, `MIMIC_PG_PORT`,
`MIMIC_PG_DBNAME`, `MIMIC_PG_USER` and `MIMIC_PG_PASSWORD` for the Python
predictor extractor. The `psql` steps use the corresponding standard libpq
environment variables (`PGHOST`, `PGPORT`, `PGDATABASE`, `PGUSER`,
`PGPASSWORD`) or the operator's secure local connection setup. Never put
credentials into scripts or commit the configured file.

Every clinical runtime output, including derived tables, logs, bootstrap
records and model objects, goes below `private_output` (default:
`results/private`, excluded from git). Runtime summaries are **not**
automatically publication-safe. `results/aggregate` is the separately screened
release snapshot and is never the destination of extraction or model fitting.

### Execution order

Run only the stages needed for a specific result. The routes below are not a
single automatic pipeline.

| Route | Order |
| --- | --- |
| MIMIC source prerequisites | Existing MIMIC derived concepts; `00_prepare_sources.sql` runs `00_source_mappings.sql`, `00_source_exclusions.sql`, `00_code_status.sql`, then `00_subject_last_record.sql` |
| MIMIC cohort and predictors | `00_extract_candidates.sql`, `00_positive_drug_check.sql`, `01_select_cohort.py`, `02_build_inputs.py` |
| MIMIC clinical results | `04_prepare_clinical.py`, then `05_cumulative_incidence.R`, `06_support_strata.R`, `07_conditional_risk.R` |
| MIMIC associations and burden | `03_association_models.R`, `08_sensitivity_models.R`; `09_prepare_burden_query.py`, execute its private SQL, then `10_restart_burden.R` |
| SICdb sources | `00_extract_source_store.py prep`, then `signals`, then `labs`; `01_select_cohort.py`; `02_build_inputs.py` |
| SICdb clinical results | `04_descriptive_tables.R`; `05_cumulative_incidence.R` with each of `S2` (overall/dose strata), `S3` (support strata), `S4` (conditional risks); `08_restart_burden.py` |
| SICdb associations | `03_association_model.R`, `06_sofa_same_subset.R`, `07_composite_model.R` |
| Amsterdam sources and cohort | `00_drug_aliases.py`; `extract_vaso_source.py`, `extract_antibiotic_process_source.py`, `extract_numeric_eligibility_source.py`, `build_technical_episodes.py`, `build_clinical_proxy_preoutcome.py`, `adjudicate_first_cohort_outcomes.py` |
| Amsterdam estimates | `03_paired_association.R`, `04_descriptive_additions.py`, `05_conditional_risk.R`, `06_restart_burden.py` |
| Common specification | `01_paired_models.R` with `mimic` and `sicdb`; Amsterdam's paired model is in its own route |
| Drug-scope comparison | `02_prepare_drug_query.py`, execute its private SQL, then `03_drug_scope.py mimic` and `03_drug_scope.py sicdb` |
| Risk model | `00_prepare_design.R`, `01_fit_ridge.R`, `02_internal_validation.R`; `03_temporal_validation.R`, `04_temporal_updates.R`; `05_sicdb_external.R`, `06_crrt_interface.R`; `07_imputation_traces.R` reads saved imputations |
| Released figures and table copies | `Rscript code/06_tables_figures/01_public_figures.R`; `python code/06_tables_figures/02_export_tables.py` |

Within each route the files are in the corresponding numbered `code/`
subdirectory. For example:

```text
psql -X -v ON_ERROR_STOP=1 -f code/01_mimic/00_prepare_sources.sql
python code/02_sicdb/00_extract_source_store.py prep
Rscript code/04_common_specification/01_paired_models.R mimic
```

For MIMIC, run the source-preparation command before
`00_extract_candidates.sql`. It creates or recreates only study-derived tables
in `project_d2v2`, `project_step4` and `project_d4`; it does not modify the MIMIC
source tables. These schemas should be dedicated to this analysis. The
candidate export is a separate read-only step.

MIMIC extraction uses `psql` output variables `candidates_file`,
`code_status_file`, `subject_last_record_file` and `positive_drug_file`.
Export these as `candidates.csv`, `code_status.csv`, `subject_last_record.csv`
and `positive_drug_rounded_zero_nee.csv` inside the configured
`mimic_candidate_exports` directory. Definitions tables are created in a local
analysis schema; source tables are read only. `00_subject_last_record.sql`
reproduces the original definition: `MAX(dischtime)` from
`mimiciv_hosp.admissions`, grouped by `subject_id`. The exported
`last_documented_dischtime` is therefore the patient's latest documented
hospital discharge, not the latest record across all clinical source tables.

Amsterdam numeric extraction expects the original inner `numericitems.zip`,
not the outer download bundle. Source/schema guards retain the version-specific
checks in the clinical implementation. The raw directory also contains
`admissions.csv`, `drugitems.csv`, `processitems.csv` and `listitems.csv` as
required by the source scripts.

## Code and manuscript map

Display CSVs mirror the manuscript table structure after disclosure protection.
Underlying numerical summaries are provided separately where safe.

| Manuscript item | Producing analysis route | Released summary |
| --- | --- | --- |
| Table 1 | Database-specific observed inputs and descriptive summaries | `Table_1_display.csv` |
| Figure 1 | Database-specific cohort selection and outcomes | `Figure_1_structure.csv` |
| Figure 2 | Clinical incidence and conditional-risk routes | `Figure_2_fixed_times.csv` |
| Figure 3 | MIMIC/SICdb support strata and conditional risks | `Figure_3_fixed_times.csv` |
| Figure 4 | MIMIC/SICdb primary associations; Amsterdam paired model; common specification | `Figure_4_source_data.csv` |
| Table S1 | Database definitions documented in the manuscript | `Table_S1_display.csv` |
| Section S2 | Analysis sequence in the manuscript supplement | Narrative only |
| Tables S3A-F | MIMIC sensitivities; SICdb subset/composite; conditional and censored incidence | `Table_S3A_display.csv` through `Table_S3F_display.csv`; association summaries |
| Tables S4A-F | Risk-model fitting and internal, temporal and external evaluation | `Table_S4A_display.csv` through `Table_S4F_display.csv`; ridge/performance/instability summaries |
| Tables S5A-H | Amsterdam clinical, burden and paired/common-specification comparisons | `Table_S5A_display.csv` through `Table_S5H_display.csv` |
| Tables S6A-B | MIMIC/SICdb baseline descriptions | `Table_S6A_display.csv`, `Table_S6B_display.csv` |
| Table S8 | Outcome-specific descriptive summaries | `Table_S8_display.csv` |
| Supplementary diagnostic figures | Risk-model calibration/stability routines and saved-imputation traces | Generated privately from local artifacts; patient-resolution inputs are not distributed |
| Supplementary Figure S6 | MIMIC/SICdb peak-dose clinical incidence | `Figure_S6_fixed_times.csv` |

`code/06_tables_figures/` draws repository derivatives directly from these
released CSVs, without fitting a model or estimating an outcome. The full
manuscript plots require the authorized local runtime outputs.

## Definition essentials

The cessation time is t0. Qualification is assessed four hours later; follow-up
ends at the original t0 plus 72 hours or the applicable observation boundary.
The primary clinical endpoint is same-care-unit resumption of the included
vasopressors. Death before restart and normal care-period ending are competing
events in the main MIMIC/SICdb cumulative-incidence presentation.
Undetermined endpoints remain separate from negative outcomes.

Definitions are database-specific. Amsterdam uses four mapped drugs,
grouped/interval inputs and recorded support processes; its fixed-window
proportions are not standard competing-risk curves. See Supplementary Table S1
for the complete three-cohort definition comparison. Cohort-specific
associations and individual-risk model transport are separate questions.

Some analyses were specified after earlier results were known; their sequence
is documented in Supplementary Section S2 of the manuscript.

## Licenses and citation

Code is MIT licensed (`LICENSE`). Released aggregate results and repository
figure derivatives are CC BY 4.0 (`LICENSE-results`). Database access and
reuse remain governed by each provider's terms. Use `CITATION.cff`; the
archival DOI and version metadata will be added when available.
