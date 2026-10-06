source("code/00_setup/paths.R")
options(stringsAsFactors=FALSE, warn=1)
suppressPackageStartupMessages({library(data.table);library(digest);library(Hmisc)})
root <- private_output("risk_model")
out <- root
input <- release_path("mimic_model_input")
source("code/05_risk_model/model_helpers.R")
dt <- fread(input)
sap_assert_main_data(dt)
allowed_missing <- c(SAP_IMPUTE_PMM,SAP_IMPUTE_LOGREG)
stopifnot(!anyNA(dt[,setdiff(SAP_PREDICTORS,allowed_missing),with=FALSE]))
stopifnot(all(vapply(dt[,SAP_PREDICTORS,with=FALSE],function(x)all(is.finite(x)|is.na(x)),logical(1))))
old_knots <- release_path("frozen_knots")
new_knots <- file.path(out,"sap3_rcs_knots_v1_3.csv")
if(file.exists(new_knots)) stopifnot(identical(digest(file=old_knots,algo="sha256"),digest(file=new_knots,algo="sha256"))) else stopifnot(file.copy(old_knots,new_knots))
knots <- sap_load_knots(new_knots)
folds <- sap_make_patient_folds(dt,seed=2026072803L)
row_fold <- sap_fold_id_for_rows(dt,folds)
fold_qc <- data.table(subject=dt$subject_id,y=dt[[SAP_OUTCOME]],fold=row_fold)[,.(landmarks=.N,patients=uniqueN(subject),events=sum(y)),by=fold][order(fold)]
stopifnot(nrow(fold_qc)==10, sum(fold_qc$landmarks)==1958, sum(fold_qc$patients)==1856,sum(fold_qc$events)==680)
# Median/mode fill is only a structural design probe, never model training.
probe <- as.data.frame(copy(dt))
for(v in allowed_missing) {
  x <- probe[[v]]
  fill <- if(v %in% SAP_BINARY_NAMES) as.integer(mean(x,na.rm=TRUE)>=0.5) else median(x,na.rm=TRUE)
  x[is.na(x)] <- fill; probe[[v]] <- x
}
x <- sap_build_design(probe,knots)
rank <- qr(cbind(1,x))$rank
condition <- kappa(scale(x),exact=TRUE)
stopifnot(ncol(x)==27,rank==28,all(apply(x,2,sd)>0))
missing <- data.table(predictor=SAP_PREDICTORS,missing_n=vapply(dt[,SAP_PREDICTORS,with=FALSE],function(x)sum(is.na(x)),integer(1)))
missing[,missing_percent:=100*missing_n/nrow(dt)]
stopifnot(all(missing$missing_percent<=40))
fwrite(folds,file.path(out,"sap3_main_patient_fold_map_v1_3.csv"))
fwrite(fold_qc,file.path(out,"sap3_main_fold_qc_v1_3.csv"))
fwrite(missing,file.path(out,"sap3_input_missingness_v1_3.csv"))
snapshot <- data.table(item=c("input_file","input_sha256","landmarks","patients","primary_events","raw_predictors","design_columns_nonintercept","matrix_rank_with_intercept","scaled_condition_number","frozen_knots_sha256","mice_m","mice_maxit","mice_pmm_donors","cv_folds","lambda_grid_points","bootstrap_success_target"),
 value=as.character(c(input,digest(file=input,algo="sha256"),1958,1856,680,13,27,rank,condition,digest(file=new_knots,algo="sha256"),SAP_M,SAP_MAXIT,SAP_DONORS,SAP_K_FOLDS,length(SAP_LAMBDA_GRID),1000)))
fwrite(snapshot,file.path(out,"sap3_input_snapshot_v1_3.csv"))
fwrite(data.table(lambda=SAP_LAMBDA_GRID),file.path(out,"sap3_lambda_grid_v1_3.csv"))
cat("PREFLIGHT_PASS: 1958 events / 1856 patients / 680 restarts; 27 design columns; original knots byte-identical.\n")
print(fold_qc)
print(missing)
cat("Scaled design condition number (structural probe):",condition,"\n")
