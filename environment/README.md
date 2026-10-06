# Environment

The recorded analysis environment is R 4.6.0. `sessionInfo.txt` contains the
available original analysis-session records, not an environment lock recreated
after the study. No original `renv.lock` was available; none has been invented.

The R code uses data.table, digest, jsonlite, mice, glmnet, Hmisc, splines,
survival, lpSolve, ggplot2 and patchwork. See the session records for recorded
versions. Graphical scripts require a working Cairo PDF device and Arial or a
compatible sans-serif font.

Python 3.11 or later is required by source extraction helpers. Install the
direct dependencies in `requirements.txt` in an isolated environment. This is
a dependency list, not a historical version lock. PostgreSQL with MIMIC-IV
derived concepts and a `psql` client is required for MIMIC extraction. SQLite is
included with Python. Packages and source databases are not bundled.

Only syntax and packaging checks have been performed on the portable copies;
the analyses have not been rerun in a clean environment.
