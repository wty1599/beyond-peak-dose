# Supplement only. Original database-specific outcomes, fixed samples and raw inputs.
source('code/00_setup/paths.R')
here <- private_output('common_specification')
stage <- file.path(here,'private'); dir.create(stage,recursive=TRUE,showWarnings=FALSE)
db <- tail(commandArgs(TRUE),1L); stopifnot(db %in% c('mimic','sicdb'))
prefix <- file.path(here,paste0(db,'_common_'))
outputs <- paste0(prefix,c('contrasts.csv','paired_changes.csv','bootstrap_replicates.csv',
                         'coefficients.csv','failures.csv','diagnostics.json'))
if(any(file.exists(outputs)))stop('refuse_output_overwrite')
stopifnot(requireNamespace('digest',quietly=TRUE),requireNamespace('jsonlite',quietly=TRUE))
sha <- function(p)tolower(digest::digest(file=p,algo='sha256'))
started <- Sys.time()
if(db=='mimic') {
  source <- release_path('mimic_observed_input')
  raw <- read.csv(source,stringsAsFactors=FALSE,check.names=FALSE)
  stopifnot(nrow(raw)==1958L,sum(raw$outcome_restart_72h)==680L,
            length(unique(raw$subject_id))==1856L)
  renal <- raw$predictor_rrt_t0
  expected_n <- 1958L; expected_y <- 680L; seed <- 20261001L
} else {
  source <- release_path('sicdb_clinical_input')
  raw <- read.csv(source,stringsAsFactors=FALSE,check.names=FALSE)
  renal <- raw$predictor_rrt_t0_recorded
  expected_n <- 1553L; expected_y <- 517L; seed <- 20261002L
}
source_hash <- sha(source)
stopifnot(all(raw$predictor_mechanical_ventilation_t0 %in% c(0,1,NA)),
          all(renal %in% c(0,1,NA)),all(raw$predictor_female %in% c(0,1,NA)))
data <- data.frame(patientid=raw$subject_id,eventid=raw$stay_id,
  y=if(db=='mimic')raw$outcome_restart_72h else raw$restart_binary,
  vent=raw$predictor_mechanical_ventilation_t0,renal=renal,
  age=cut(raw$predictor_age,breaks=c(18,40,50,60,70,80,Inf),right=FALSE,
          labels=c('18-39','40-49','50-59','60-69','70-79','80+')),
  sex=factor(ifelse(raw$predictor_female==1,'Vrouw',ifelse(raw$predictor_female==0,'Man',NA)),
             levels=c('Vrouw','Man')),
  log2duration=log2(raw$predictor_positive_nee_hours),
  peak_band=cut(raw$predictor_episode_peak_nee,breaks=c(-Inf,.2,.5,Inf),right=FALSE,
                labels=c('<0.20','0.20-<0.50','>=0.50')))
complete <- complete.cases(data)&is.finite(data$log2duration)
attrition <- list(source_events=nrow(data),original_restart_unknown=sum(is.na(data$y)),
                 any_required_input_or_outcome_missing=sum(!complete),
                 ventilation_unknown=sum(is.na(data$vent)),renal_unknown=sum(is.na(data$renal)),
                 age_unknown=sum(is.na(data$age)),sex_unknown=sum(is.na(data$sex)),
                 peak_band_unknown=sum(is.na(data$peak_band)),duration_invalid=sum(!is.finite(data$log2duration)))
data <- data[complete,,drop=FALSE]
stopifnot(nrow(data)==expected_n,sum(data$y)==expected_y,!anyDuplicated(data$eventid),
          all(data$y %in% 0:1),all(table(data$age)>0),all(table(data$sex)>0),all(table(data$peak_band)>0))
if(db=='sicdb') {
  s30source <- release_path('sicdb_common_identity_input')
  s30 <- read.csv(s30source,stringsAsFactors=FALSE)
  stopifnot(nrow(s30)==1553L,setequal(data$eventid,s30$stay_id),
            all(raw$predictor_rrt_t0_recorded[match(s30$CaseID,raw$CaseID)]==s30$recorded_crrt_active_pre6h))
}
forms <- list(without_dose_band=y~vent+renal+age+sex+log2duration,
              with_dose_band=y~vent+renal+age+sex+log2duration+peak_band)
fit_checked <- function(dat,form) {
  warn <- character()
  fit <- withCallingHandlers(glm(form,data=dat,family=binomial(),control=glm.control(maxit=50L)),
    warning=function(w){warn <<- c(warn,conditionMessage(w));invokeRestart('muffleWarning')})
  cf <- coef(fit); p <- fitted(fit)
  if(length(warn)||!fit$converged||fit$rank!=length(cf)||any(!is.finite(cf))||
     any(!is.finite(p))||any(p<=1e-7|p>=1-1e-7))
    stop(paste(c('fit_gate_failure',warn),collapse=';'))
  fit
}
estimate <- function(dat,fit) {
  out <- numeric()
  for(term in c('vent','renal')) {
    zero <- one <- dat; zero[[term]]<-0L; one[[term]]<-1L
    r0 <- mean(predict(fit,newdata=zero,type='response'))
    r1 <- mean(predict(fit,newdata=one,type='response'))
    out[paste0(term,'_risk_no')]<-r0;out[paste0(term,'_risk_yes')]<-r1
    out[paste0(term,'_rd_pp')]<-100*(r1-r0);out[paste0(term,'_or')]<-exp(coef(fit)[[term]])
  }
  stopifnot(all(is.finite(out)));out
}
fits <- lapply(forms,function(f)fit_checked(data,f))
stopifnot(length(coef(fits$without_dose_band))==10L,length(coef(fits$with_dose_band))==12L)
points <- lapply(fits,function(f)estimate(data,f))
clusters <- split(seq_len(nrow(data)),data$patientid);ids <- names(clusters)
B <- 1000L;set.seed(seed)
boot <- array(NA_real_,dim=c(B,8L,2L),dimnames=list(NULL,names(points[[1]]),names(forms)))
coeff <- lapply(fits,function(f)matrix(NA_real_,nrow=B,ncol=length(coef(f)),dimnames=list(NULL,names(coef(f)))))
membership <- vector('list',B);failures <- list()
for(b in seq_len(B)) {
  drawn <- sample(ids,length(ids),replace=TRUE);index <- unlist(clusters[drawn],use.names=FALSE)
  membership[[b]] <- index;sample_data <- data[index,,drop=FALSE]
  for(m in names(forms)) {
    tryCatch({f<-fit_checked(sample_data,forms[[m]])
      boot[b,,m]<-estimate(sample_data,f);coeff[[m]][b,]<-coef(f)},
      error=function(e){failures[[length(failures)+1L]] <<-
        data.frame(draw=b,model=m,reason=conditionMessage(e))})
  }
  if(b%%100L==0L)cat(db,': paired draws completed ',b,'/',B,'\n',sep='')
}
valid <- vapply(seq_len(B),function(b)all(is.finite(boot[b,,])),logical(1))
can_ci <- sum(valid)>=990L
ci <- function(x)if(can_ci)as.numeric(quantile(x[valid],c(.025,.975),type=7))else c(NA_real_,NA_real_)
rows <- list();k<-0L
for(m in names(forms))for(metric in names(points[[m]])) {
  lim<-ci(boot[,metric,m]);k<-k+1L
  rows[[k]]<-data.frame(database=if(db=='mimic')'MIMIC-IV' else 'SICdb',model=m,metric=metric,
    estimate=unname(points[[m]][metric]),ci_low=lim[1],ci_high=lim[2],
    n_events=nrow(data),n_patients=length(ids),n_restarts=sum(data$y),
    paired_bootstrap_success=sum(valid),paired_bootstrap_failure=sum(!valid))
}
write.csv(do.call(rbind,rows),outputs[1],row.names=FALSE,na='')
change_rows <- lapply(c('vent','renal'),function(term){metric<-paste0(term,'_rd_pp')
  delta<-boot[,metric,'with_dose_band']-boot[,metric,'without_dose_band'];lim<-ci(delta)
  data.frame(database=if(db=='mimic')'MIMIC-IV' else 'SICdb',indicator=term,
    contrast='with_minus_without_dose_band',rd_change_pp=unname(points$with_dose_band[metric]-points$without_dose_band[metric]),
    ci_low=lim[1],ci_high=lim[2],n_events=nrow(data),n_patients=length(ids),n_restarts=sum(data$y),
    paired_bootstrap_success=sum(valid),paired_bootstrap_failure=sum(!valid))})
write.csv(do.call(rbind,change_rows),outputs[2],row.names=FALSE,na='')
replicates <- data.frame(draw=seq_len(B),jointly_valid=valid)
for(m in names(forms))for(metric in names(points[[m]]))replicates[[paste(m,metric,sep='__')]]<-boot[,metric,m]
write.csv(replicates,outputs[3],row.names=FALSE,na='')
coef_rows <- list();k<-0L
for(m in names(forms))for(term in names(coef(fits[[m]]))){lim<-ci(coeff[[m]][,term]);k<-k+1L
  coef_rows[[k]]<-data.frame(model=m,term=term,log_odds_coefficient=unname(coef(fits[[m]])[term]),ci_low=lim[1],ci_high=lim[2])}
write.csv(do.call(rbind,coef_rows),outputs[4],row.names=FALSE,na='')
fail_table <- if(length(failures))do.call(rbind,failures)else data.frame(draw=integer(),model=character(),reason=character())
write.csv(fail_table,outputs[5],row.names=FALSE)
saveRDS(list(data=data,membership=membership,forms=forms,seed=seed,points=points,boot=boot,
             coefficient_boot=coeff,original_source=source),file.path(stage,paste0(db,'_common_reproduction.rds')),compress='gzip')
diagnostics <- lapply(fits,function(f)list(converged=f$converged,rank=f$rank,parameter_count=length(coef(f)),
  fitted_probability_range=range(fitted(f)),largest_10_cook_distances=unname(head(sort(cooks.distance(f),decreasing=TRUE),10))))
stopifnot(sha(source)==source_hash)
audit <- list(started=as.character(started),finished=as.character(Sys.time()),
  elapsed_seconds=as.numeric(difftime(Sys.time(),started,units='secs')),R_version=R.version.string,
  package_versions=list(stats=as.character(packageVersion('stats')),digest=as.character(packageVersion('digest')),jsonlite=as.character(packageVersion('jsonlite'))),
  seed=seed,bootstrap_draws=1000L,events=nrow(data),patients=length(ids),restarts=sum(data$y),
  jointly_successful=sum(valid),jointly_failed=sum(!valid),interpretation_gate_pass=can_ci,
  same_rows_both_models=TRUE,standardisation='each bootstrap draw sampled events',percentile_type=7,
  original_outcome_preserved=TRUE,no_winsorisation_or_spline_reselection=TRUE,source_sha256=source_hash,
  source_unchanged=TRUE,attrition=attrition,age_counts=as.list(table(data$age)),sex_counts=as.list(table(data$sex)),
  peak_band_counts=as.list(table(data$peak_band)),diagnostics=diagnostics,no_individual_predictions_exported=TRUE)
jsonlite::write_json(audit,outputs[6],auto_unbox=TRUE,pretty=TRUE,digits=NA)
cat('Completed ',db,': ',sum(valid),' jointly successful draws.\n',sep='')
if(!can_ci)stop('More_than_one_percent_failed; report only, no interpretation')
