options(stringsAsFactors=FALSE,warn=1)
suppressPackageStartupMessages({library(data.table);library(mice);library(Hmisc);library(digest);library(jsonlite)})
source('code/00_setup/paths.R')
root<-private_output('sicdb_external');a<-file.path(root,'aggregate');w<-file.path(root,'private')
dir.create(a,recursive=TRUE,showWarnings=FALSE);dir.create(w,recursive=TRUE,showWarnings=FALSE)
paths<-c(train=release_path('mimic_model_input'),artifact=file.path(private_output('risk_model'),'sap3_main_model_artifact_v1_3.rds'),helper='code/05_risk_model/model_helpers.R',external=release_path('sicdb_scorable_input'),winsor=release_path('frozen_winsor'),knots=release_path('frozen_knots'))
hashes<-vapply(paths,function(p)digest(file=p,algo='sha256'),character(1))
source(paths['helper'])
train <- fread(paths['train']);sap_assert_main_data(train)
artifact <- readRDS(paths['artifact'])
stopifnot(artifact$protocol_version=='original_method_minimal_correction_20260920',artifact$input_sha256==hashes['train'],
  identical(artifact$design_columns,sap_design_column_names),identical(dim(artifact$coefficient_matrix),c(28L,20L)),
  isTRUE(all.equal(artifact$knots,sap_load_knots(paths['knots']),check.attributes=TRUE)))
raw <- fread(paths['external']);external <- raw[,c(SAP_ID_NAMES,SAP_OUTCOME,SAP_PREDICTORS),with=FALSE]
y <- external[[SAP_OUTCOME]];n <- nrow(external);np <- uniqueN(external$subject_id)
stopifnot(n==789L,np==772L,sum(y)==277L,all(y%in%c(0,1)),uniqueN(external$model_row_id)==n)
required <- setdiff(SAP_PREDICTORS,c(SAP_IMPUTE_PMM,SAP_IMPUTE_LOGREG))
stopifnot(!anyNA(external[,..required]))
win <- fread(paths['winsor'])
for(i in seq_len(nrow(win))){v <- win$predictor[i];set(external,j=v,value=pmin(pmax(external[[v]],win$lower_p005[i]),win$upper_p995[i]))}
missing <- data.table(predictor=SAP_PREDICTORS,n=n,missing=vapply(external[,..SAP_PREDICTORS],function(z)sum(is.na(z)),integer(1)))
stopifnot(sum(missing$missing)==20L,missing[predictor=='predictor_urine_output_24h_ml',missing]==20L)
fwrite(missing,file.path(a,'scorable_missingness_before_test_imputation.csv'))
fwrite(data.table(role=names(paths),path=unname(paths),sha256=unname(hashes)),file.path(a,'input_source_manifest.csv'))
write_json(list(recorded_at=as.character(Sys.time()),model=artifact$protocol_version,n=n,patients=np,events=sum(y),training_n=nrow(train),
  phase='BEFORE_EXTERNAL_PREDICTION',model_refit=FALSE,cv=FALSE,seed_mice=202609091L,seed_bootstrap=202609092L,
  m=20L,maxit=10L,external_outcomes_in_imputer=FALSE,neutral_y=mean(train[[SAP_OUTCOME]]),
  coefficient_artifact_sha256=hashes['artifact'],external_input_sha256=hashes['external'],
  strict_external_validation=FALSE,recalibration_updates=FALSE),file.path(a,'prediction_start_record.json'),pretty=TRUE,auto_unbox=TRUE)
combined <- rbindlist(list(train[,c(SAP_OUTCOME,SAP_PREDICTORS),with=FALSE],external[,c(SAP_OUTCOME,SAP_PREDICTORS),with=FALSE]))
combined[,(SAP_OUTCOME):=as.numeric(get(SAP_OUTCOME))]
test <- seq.int(nrow(train)+1L,nrow(combined));combined[test,(SAP_OUTCOME):=mean(train[[SAP_OUTCOME]])]
cat('EXTERNAL_TEST_IMPUTATION_ONLY n=789; FROZEN_MODEL n=1958\n');flush.console()
mi <- sap_run_mice(combined,seed=202609091L,ignore=seq_len(nrow(combined))%in%test,neutral_outcome_rows=test,m=20L,maxit=10L)
completed <- lapply(mi$completed,function(z)z[test,,drop=FALSE])
pr <- sap_predict_completed_from_coefficients(completed,artifact$coefficient_matrix,artifact$knots)
p <- pr$pooled_probability
stopifnot(length(p)==789L,ncol(pr$predictions)==20L,all(is.finite(p)),all(p>0&p<1))
pred <- external[,c(SAP_ID_NAMES,SAP_OUTCOME),with=FALSE];pred[,probability:=p]
fwrite(pred,file.path(w,'external_predictions.csv'))
saveRDS(list(predictions=pr,completed=completed,model=artifact$protocol_version,source_sha256=hashes,mice_methods=mi$method,predictor_matrix=mi$predictor_matrix),file.path(w,'external_prediction_artifact.rds'))
saveRDS(mi$mids,file.path(w,'external_test_imputation_mids.rds'))
fwrite(sap_mice_log_table(mi,'SICdb_current1958_test_ignore'),file.path(a,'test_mice_logged_events.csv'))
metrics <- function(yy,pp){
  lp <- qlogis(pmin(pmax(pp,1e-8),1-1e-8))
  off <- glm(yy~1,offset=lp,family=binomial());joint <- glm(yy~lp,family=binomial())
  stopifnot(off$converged,joint$converged)
  c(c_statistic=sap_auc(yy,pp),calibration_intercept_offset=unname(coef(off)[1]),
    calibration_slope_joint=unname(coef(joint)[2]),calibration_intercept_joint=unname(coef(joint)[1]),
    brier=mean((yy-pp)^2),observed_risk=mean(yy),mean_prediction=mean(pp))
}
original <- metrics(y,p)
groups <- split(seq_len(n),external$subject_id);set.seed(202609092L)
bs <- vector('list',1000L);failed <- list();success <- 0L;attempt <- 0L
while(success<1000L && attempt<1500L){
  attempt <- attempt+1L;ii <- unlist(groups[sample.int(length(groups),length(groups),replace=TRUE)],use.names=FALSE)
  v <- tryCatch({if(length(unique(y[ii]))!=2L)stop('one outcome class');z <- metrics(y[ii],p[ii]);if(any(!is.finite(z)))stop('nonfinite');z},error=function(e)e)
  if(inherits(v,'error')) failed[[length(failed)+1L]] <- data.table(attempt=attempt,error=conditionMessage(v)) else {success <- success+1L;bs[[success]] <- c(attempt=attempt,v)}
}
stopifnot(success==1000L)
bt <- rbindlist(lapply(bs,as.list));fwrite(bt,file.path(a,'external_fixed_prediction_bootstrap.csv'))
ff <- if(length(failed)) rbindlist(failed) else data.table(attempt=integer(),error=character())
fwrite(ff,file.path(a,'external_bootstrap_failures.csv'))
perf <- rbindlist(lapply(names(original),function(k)data.table(model=artifact$protocol_version,metric=k,estimate=unname(original[k]),lower95=quantile(bt[[k]],.025),upper95=quantile(bt[[k]],.975),n=n,patients=np,events=sum(y),interval='Patient-cluster percentile; fixed external predictions')))
fwrite(perf,file.path(a,'external_performance_current_model.csv'))
stopifnot(identical(hashes,vapply(paths,function(p)digest(file=p,algo='sha256'),character(1))))
write_json(list(status='PASS_ADAPTED_EXTERNAL_EVALUATION',model=artifact$protocol_version,n=n,patients=np,events=sum(y),
  source_files_unchanged=TRUE,protected_identity_files=31L,model_refits=0L,lambda_reselection=0L,cv_runs=0L,
  development_mids_reruns=0L,external_test_imputation_runs=1L,external_ignore_rows=length(test),external_true_outcomes_in_imputer=FALSE,
  probability_aggregation='mean_of_20_probabilities',imputation_logged_events=mi$logged_count,
  performance_bootstrap_success=success,performance_bootstrap_attempts=attempt,performance_bootstrap_failures=nrow(ff),
  recalibration_updates=0L,strict_external_validation=FALSE,CRRT_original_case_level_proxy_limitation_retained=TRUE),file.path(a,'external_validation_qc.json'),pretty=TRUE,auto_unbox=TRUE)
capture.output(sessionInfo(),file=file.path(a,'R_session_info.txt'))
print(perf)
cat('CURRENT1958_MODEL_EXTERNAL_EVALUATION_COMPLETE\n')
