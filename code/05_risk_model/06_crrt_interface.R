options(stringsAsFactors=FALSE,digits=17,warn=1)
suppressPackageStartupMessages({library(data.table);library(digest);library(jsonlite)})
source('code/00_setup/paths.R')
out<-private_output('sicdb_external_crrt');start<-Sys.time()
ext<-private_output('sicdb_external')
paths<-c(external=release_path('sicdb_scorable_input'),s30=release_path('sicdb_association_input'),baseline=file.path(ext,'aggregate/external_performance_current_model.csv'),boots=file.path(ext,'aggregate/external_fixed_prediction_bootstrap.csv'),predictions=file.path(ext,'private/external_predictions.csv'),model=file.path(private_output('risk_model'),'sap3_main_model_artifact_v1_3.rds'))
hash <- vapply(paths,function(p)digest(file=p,algo='sha256'),character(1))
d <- fread(paths['external']);s <- fread(paths['s30'])
stopifnot(nrow(d)==789,uniqueN(d$PatientID)==772,uniqueN(d$CaseID)==789)
rrt <- d$recorded_crrt_active_pre6h
qc <- data.table(check=c('scorable_events','recorded_crrt_missing','recorded_crrt_not_binary','legacy_vs_recorded_changed','matches_S30_cases'),
 value=c(nrow(d),sum(is.na(rrt)),sum(!rrt%in%0:1),sum(d$predictor_rrt_t0!=rrt,na.rm=TRUE),sum(d$CaseID%in%s$CaseID)))
fwrite(qc,file.path(out,'interface_checks.csv'))
if(anyNA(rrt)||any(!rrt%in%0:1)) stop('CRRT interface unclassifiable: await author')
ix <- match(d$CaseID,s$CaseID)
stopifnot(!anyNA(ix),all(d$t0==s$t0[ix]),all(rrt==s$predictor_rrt_t0[ix]))
if(any(d$predictor_rrt_t0!=rrt))stop('INPUT_CHANGED: requires same original deployment imputation, not shortcut')
# existing 20-stream completed data and predictions exact reusable results.
p <- fread(paths['predictions']);stopifnot(identical(p$model_row_id,d$model_row_id))
yy <- p$outcome_restart_72h;pp <- p$probability;lp <- qlogis(pp)
off <- glm(yy~1,offset=lp,family=binomial());joint <- glm(yy~lp,family=binomial())
stopifnot(off$converged,joint$converged)
n1 <- sum(yy==1);n0 <- sum(yy==0)
auc <- (sum(rank(pp,ties.method='average')[yy==1])-n1*(n1+1)/2)/(n1*n0)
v <- c(c_statistic=auc,calibration_intercept_offset=unname(coef(off)[1]),calibration_slope_joint=unname(coef(joint)[2]),calibration_intercept_joint=unname(coef(joint)[1]),brier=mean((yy-pp)^2),observed_risk=mean(yy),mean_prediction=mean(pp))
base <- fread(paths['baseline']);boot <- fread(paths['boots']);stopifnot(nrow(boot)==1000)
delta <- v[base$metric]-base$estimate
stopifnot(max(abs(delta))<1e-12)
for(i in seq_len(nrow(base))){ci <- quantile(boot[[base$metric[i]]],c(.025,.975));stopifnot(max(abs(ci-c(base$lower95[i],base$upper95[i])))<1e-12)}
result <- rbindlist(list(copy(base)[,interface:='legacy_case_level'],copy(base)[,interface:='recorded_crrt_pre6h']))
result[,basis:='Identical 789 inputs; reuse same 20 MI predictions and same 1000 patient-cluster draws']
fwrite(result,file.path(out,'T4B_performance_side_by_side.csv'))
fwrite(data.table(metric=names(v),independent_estimate=unname(v),source_estimate=base$estimate[match(names(v),base$metric)],absolute_difference=abs(unname(v)-base$estimate[match(names(v),base$metric)])),file.path(out,'T4B_numeric_check.csv'))
after <- vapply(paths,function(p)digest(file=p,algo='sha256'),character(1));stopifnot(identical(hash,after))
fwrite(data.table(role=names(paths),path=unname(paths),sha256_before=unname(hash),sha256_after=unname(after),unchanged=hash==after),file.path(out,'source_protection.csv'))
write_json(list(status='PASS_IDENTICAL_INTERFACE_ON_FIXED_SUBSET',n=789,patients=772,changed_inputs=0,unknown_recorded_interface=0,imputations_preserved=20,bootstrap_draws_preserved=1000,bootstrap_seed=202609092,imputation_seed=202609091,new_imputations=0,new_prediction_model_fits=0,limitation='Does not remove prior case-level selection or establish absence of all RRT'),file.path(out,'T4B_QC.json'),pretty=TRUE,auto_unbox=TRUE)
cat('Elapsed seconds:',as.numeric(difftime(Sys.time(),start,units='secs')),'\n');print(sessionInfo());
