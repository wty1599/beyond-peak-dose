options(stringsAsFactors=FALSE,warn=1)
suppressPackageStartupMessages({library(data.table);library(Hmisc);library(jsonlite)})
source('code/00_setup/paths.R')
root <- private_output('sicdb_association')
a <- root; w <- file.path(root,'private')
dir.create(w,recursive=TRUE,showWarnings=FALSE)
if(file.exists(file.path(a,'association_model_qc.json')) || file.exists(file.path(a,'association_cluster_bootstrap.csv'))) stop('Preserve completed or partial analysis')
source('code/05_risk_model/model_helpers.R')
v <- fread(release_path('sicdb_association_input'))
gate <- list(n=nrow(v), patients=uniqueN(v$PatientID), events=sum(v$outcome_restart_72h), events_at_least_400=sum(v$outcome_restart_72h)>=400, duration_IQR_at_least_30h=IQR(v$predictor_positive_nee_hours)>=30)
stopifnot(nrow(v)==gate$n,uniqueN(v$PatientID)==gate$patients,sum(v$outcome_restart_72h)==gate$events)
req <- c('predictor_age','predictor_female','predictor_episode_peak_nee','predictor_positive_nee_hours','predictor_mechanical_ventilation_t0','predictor_rrt_t0')
stopifnot(!anyNA(v[,..req]),!anyDuplicated(v$CaseID),!('predictor_sofa_24h'%in%names(v)))
knots <- sap_load_knots(release_path('frozen_knots'))
win <- fread(release_path('frozen_winsor'))
raw <- copy(v)
clipcounts <- list()
for(i in seq_len(nrow(win))){k<-win$predictor[i];if(k%in%req){
 clipcounts[[length(clipcounts)+1L]]<-data.table(variable=k,below_lower=sum(v[[k]]<win$lower_p005[i]),above_upper=sum(v[[k]]>win$upper_p995[i]),lower=win$lower_p005[i],upper=win$upper_p995[i])
 set(v,j=k,value=pmin(pmax(v[[k]],win$lower_p005[i]),win$upper_p995[i]))}}
fwrite(rbindlist(clipcounts),file.path(a,'association_frozen_winsor_application.csv'))
basis <- function(k,z){m<-Hmisc::rcspline.eval(z,knots=knots[[k]],inclx=TRUE,norm=2);colnames(m)<-paste0(k,'_rcs',1:3);m}
x <- cbind('(Intercept)'=1,basis('predictor_age',v$predictor_age),predictor_female=v$predictor_female,
 basis('predictor_episode_peak_nee',v$predictor_episode_peak_nee),basis('predictor_positive_nee_hours',v$predictor_positive_nee_hours),
 predictor_mechanical_ventilation_t0=v$predictor_mechanical_ventilation_t0,predictor_rrt_t0=v$predictor_rrt_t0)
y <- v$outcome_restart_72h
stopifnot(ncol(x)==13L,qr(x)$rank==13L,!anyNA(x),all(is.finite(x)),all(y%in%c(0,1)))
fwrite(data.table(column=colnames(x),position=seq_len(ncol(x))),file.path(a,'association_design_columns.csv'))
exposures <- c('predictor_positive_nee_hours','predictor_episode_peak_nee','predictor_mechanical_ventilation_t0','predictor_rrt_t0')
values <- list(c(24,48),c(.15,.30),c(0,1),c(0,1))
grids <- lapply(exposures[1:2],function(k)seq(quantile(raw[[k]],.05),quantile(raw[[k]],.95),length.out=40L))
clipvalue <- function(k,z){j<-match(k,win$predictor);if(!is.na(j))pmin(pmax(z,win$lower_p005[j]),win$upper_p995[j]) else z}
cols <- lapply(exposures,function(k)which(startsWith(colnames(x),k)))
newmat <- function(i,z)if(i<=2L)basis(exposures[i],clipvalue(exposures[i],z)) else matrix(z,ncol=1L)
gridbasis <- lapply(1:2,function(i)newmat(i,grids[[i]]));contrastbasis<-lapply(1:4,function(i)newmat(i,values[[i]]))
auc <- function(y,p){n1<-sum(y==1);n0<-sum(y==0);if(min(n1,n0)==0L)stop('No events or non-events');(sum(rank(p,ties.method='average')[y==1])-n1*(n1+1)/2)/(n1*n0)}
fitset <- function(xx,yy){
 jj<-list(full=seq_len(ncol(xx)),drop_duration=setdiff(seq_len(ncol(xx)),cols[[1]]),drop_peak=setdiff(seq_len(ncol(xx)),cols[[2]]))
 f<-lapply(jj,function(j){
  warnings_seen<-character()
  z<-withCallingHandlers(glm.fit(xx[,j,drop=FALSE],yy,family=binomial(),control=glm.control(maxit=100,epsilon=1e-8)),warning=function(e){warnings_seen<<-c(warnings_seen,conditionMessage(e));invokeRestart('muffleWarning')})
  if(!z$converged||z$rank!=length(j)||any(!is.finite(z$coefficients)))stop('nonconvergence/rank/coefficients')
  z$warnings_seen<-warnings_seen;z})
 cc<-vapply(f,function(z)auc(yy,z$fitted.values),numeric(1))
 dd<-f$drop_duration$deviance-f$full$deviance;dp<-f$drop_peak$deviance-f$full$deviance
 if(dd < -1e-7 || dp < -1e-7)stop('Reduced deviance below full model')
 metrics<-c(partial_deviance_duration=dd,partial_deviance_peak=dp,delta_partial_deviance=dd-dp,
 C_full_apparent=unname(cc[1]),C_drop_duration_apparent=unname(cc[2]),C_drop_peak_apparent=unname(cc[3]),
 C_loss_duration_apparent=unname(cc[1]-cc[2]),C_loss_peak_apparent=unname(cc[1]-cc[3]),delta_C_loss_duration_minus_peak_apparent=unname(cc[3]-cc[2]))
 beta<-f$full$coefficients;lp<-drop(xx%*%beta)
 risks<-function(i,nb){base<-lp-drop(xx[,cols[[i]],drop=FALSE]%*%beta[cols[[i]]]);add<-drop(nb%*%beta[cols[[i]]]);colMeans(plogis(outer(base,add,'+')))}
 gc<-unlist(lapply(1:2,function(i)risks(i,gridbasis[[i]])))
 cv<-unlist(lapply(1:4,function(i){pp<-risks(i,contrastbasis[[i]]);c(risk_low=pp[1],risk_high=pp[2],risk_difference=pp[2]-pp[1],log_OR=drop((contrastbasis[[i]][2,,drop=FALSE]-contrastbasis[[i]][1,,drop=FALSE])%*%beta[cols[[i]]]))}))
 out<-c(metrics,setNames(gc,paste0('curve_',seq_along(gc))),setNames(cv,paste0(rep(exposures,each=4),'_',rep(c('risk_low','risk_high','risk_difference','log_OR'),4))))
 if(any(!is.finite(out)))stop('Nonfinite standardization')
 list(metrics=out,fits=f,warnings=unlist(lapply(f,function(z)z$warnings_seen)))
}
pointfit <- fitset(x,y);point <- pointfit$metrics
coefout<-rbindlist(lapply(names(pointfit$fits),function(nm){z<-pointfit$fits[[nm]];data.table(model=nm,term=names(z$coefficients),estimate=z$coefficients)}))
fwrite(coefout,file.path(a,'association_new_association_coefficients_NOT_FROZEN.csv'))
designqc<-rbindlist(lapply(names(pointfit$fits),function(nm){z<-pointfit$fits[[nm]];data.table(model=nm,n=length(y),rank=z$rank,predictor_df=z$rank-1L,converged=z$converged,iterations=z$iter,deviance=z$deviance,warnings=paste(z$warnings_seen,collapse=' | '),min_probability=min(z$fitted.values),max_probability=max(z$fitted.values))}))
fwrite(designqc,file.path(a,'association_point_model_diagnostics.csv'))
ids <- split(seq_len(nrow(v)),v$PatientID);SEED<-2026091030L;set.seed(SEED)
success<-attempt<-0L;boot<-list();fail<-list();warninglog<-list()
while(success<1000L && attempt<1500L){
 attempt<-attempt+1L;ii<-unlist(ids[sample.int(length(ids),length(ids),replace=TRUE)],use.names=FALSE)
 rr<-tryCatch(fitset(x[ii,,drop=FALSE],y[ii]),error=function(e)e)
 if(inherits(rr,'error')){fail[[length(fail)+1L]]<-data.table(attempt=attempt,reason=conditionMessage(rr))
 } else {success<-success+1L;boot[[success]]<-as.data.table(c(list(attempt=attempt),as.list(rr$metrics)))
  if(length(rr$warnings))warninglog[[length(warninglog)+1L]]<-data.table(attempt=attempt,warning=unique(rr$warnings))}
 if(attempt%%100L==0L){cat('association_BOOTSTRAP success=',success,' attempts=',attempt,'\n');flush.console()}
}
bt<-rbindlist(boot);ff<-if(length(fail))rbindlist(fail) else data.table(attempt=integer(),reason=character())
wl<-if(length(warninglog))rbindlist(warninglog) else data.table(attempt=integer(),warning=character())
fwrite(bt,file.path(a,'association_cluster_bootstrap.csv'));fwrite(ff,file.path(a,'association_bootstrap_failures.csv'));fwrite(wl,file.path(a,'association_bootstrap_warnings.csv'))
if(success!=1000L)stop('complete-case bootstrap target not met; retained partial results and failures')
ci<-function(nm)c(estimate=unname(point[nm]),lower95=unname(quantile(bt[[nm]],.025)),upper95=unname(quantile(bt[[nm]],.975)))
metrics<-rbindlist(lapply(names(point)[1:9],function(nm)as.data.table(c(list(metric=nm),as.list(ci(nm))))))
fwrite(metrics,file.path(a,'association_block_comparisons.csv'))
gc<-rbindlist(lapply(1:2,function(i)rbindlist(lapply(1:40,function(j)as.data.table(c(list(variable=exposures[i],value=grids[[i]][j]),as.list(ci(paste0('curve_',(i-1)*40+j)))))))))
fwrite(gc,file.path(a,'association_standardized_curves.csv'))
ct<-rbindlist(lapply(1:4,function(i)rbindlist(lapply(c('risk_low','risk_high','risk_difference','log_OR'),function(m){vv<-ci(paste0(exposures[i],'_',m));if(m=='log_OR')vv<-exp(vv);as.data.table(c(list(variable=exposures[i],reference=values[[i]][1],comparison=values[[i]][2],measure=if(m=='log_OR')'OR' else m),as.list(vv)))}))))
fwrite(ct,file.path(a,'association_absolute_risks_and_OR.csv'))
saveRDS(list(design_columns=colnames(x),knots=knots,exposures=exposures,values=values,standardization='empirical complete-case complete required-input cohort; bootstrap sample empirical distribution',point=point,coefficients=pointfit$fits$full$coefficients),file.path(w,'association_association_design_NOT_FROZEN_MODEL.rds'))
qc<-list(n=nrow(v),patients=uniqueN(v$PatientID),events=sum(y),predictor_df=12L,design_columns=13L,design_rank=qr(x)$rank,
 bootstrap_success=success,bootstrap_attempts=attempt,bootstrap_failures=nrow(ff),bootstrap_warning_rows=nrow(wl),seed=SEED,
 full_drop_models_same_rows_and_same_patient_resamples=TRUE,standardization='case-weighted empirical complete-case cohort; resampled cohort per bootstrap',
 event_gate_at_least_400=gate$events_at_least_400,duration_IQR_gate_at_least_30h=gate$duration_IQR_at_least_30h,
 duration_superiority_interpretation_allowed=gate$events_at_least_400 && gate$duration_IQR_at_least_30h,
 SOFA_in_design=FALSE,new_imputations=0L,frozen_prediction_model_modified=FALSE,old_models_refitted=FALSE,
 C_is_apparent_not_validation=TRUE,cancelled_AB_used=FALSE,fit_spec='independent unpenalized logistic, fixed MIMIC winsor and RCS4 knots (3 df per spline)',
 post_hoc=TRUE)
write_json(qc,file.path(a,'association_model_qc.json'),pretty=TRUE,auto_unbox=TRUE)
capture.output(sessionInfo(),file=file.path(a,'association_R_session_info.txt'))
print(metrics);print(ct);print(qc)
