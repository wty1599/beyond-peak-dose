source('code/00_setup/paths.R')
here <- private_output('amsterdam_association')
stage <- file.path(here, 'private')
dir.create(stage,recursive=TRUE,showWarnings=FALSE)
cohort <- release_path('amsterdam_cohort_input')
admissions <- release_path('amsterdam_admissions')
outputs <- file.path(here, c('contrasts.csv','paired_changes.csv','bootstrap_replicates.csv',
  'coefficients.csv','failures.csv','diagnostics.json'))
if (any(file.exists(outputs))) stop('Refusing to overwrite existing outputs')
started <- Sys.time()
d <- read.csv(gzfile(cohort),stringsAsFactors=FALSE,check.names=FALSE)
a <- read.csv(admissions,fileEncoding="CP1252",stringsAsFactors=FALSE)
stopifnot(nrow(d)==779L,length(unique(d$admissionid))==779L,
  length(unique(d$patientid))==764L,nrow(a)==23106L,!anyDuplicated(a$admissionid))
j <- match(d$admissionid,a$admissionid); stopifnot(!anyNA(j))
levels_age <- c("18-39","40-49","50-59","60-69","70-79","80+")
d$age <- factor(a$agegroup[j],levels=levels_age)
d$sex <- factor(ifelse(a$gender[j]=="",NA_character_,a$gender[j]),levels=c("Vrouw","Man"))
d$y <- ifelse(d$h72_restart=="YES",1L,ifelse(d$h72_restart=="NO",0L,NA_integer_))
stopifnot(sum(d$y==1L,na.rm=TRUE)==206L,sum(d$y==0L,na.rm=TRUE)==560L,sum(is.na(d$y))==13L)
allowed <- c("RECORDED_ACTIVE_PROCESS","NO_MAPPED_PROCESS_RECORDED")
stopifnot(all(d$recorded_ventilation_process_at_t0 %in% allowed),
  all(d$recorded_cvvh_process_at_t0 %in% allowed))
d$vent <- as.integer(d$recorded_ventilation_process_at_t0==allowed[1])
d$cvvh <- as.integer(d$recorded_cvvh_process_at_t0==allowed[1])
hours <- as.numeric(d$positive_hours_segmentwise)
stopifnot(all(is.finite(hours)),all(hours>0)); d$log2duration <- log2(hours)
lo <- as.numeric(d$four_drug_nee_peak_lower); hi <- as.numeric(d$four_drug_nee_peak_upper)
stopifnot(!any(!is.na(lo)&!is.na(hi)&lo>hi))
band <- rep(NA_character_,nrow(d))
band[!is.na(hi)&hi<.20] <- "<0.20"
band[!is.na(lo)&!is.na(hi)&lo>=.20&hi<.50] <- "0.20-<0.50"
band[!is.na(lo)&lo>=.50] <- ">=0.50"
d$peak_band <- factor(band,levels=c("<0.20","0.20-<0.50",">=0.50"))
primary <- d[!is.na(d$y)&!is.na(d$age)&!is.na(d$sex)&is.finite(d$log2duration),,drop=FALSE]
data <- primary[!is.na(primary$peak_band),,drop=FALSE]
stopifnot(nrow(primary)==756L,nrow(data)==634L,length(unique(data$patientid))==624L,sum(data$y)==181L)
forms <- list(without_dose_band=y~vent+cvvh+age+sex+log2duration,
              with_dose_band=y~vent+cvvh+age+sex+log2duration+peak_band)
fit_checked <- function(dat,form) {
  warn <- character()
  fit <- withCallingHandlers(glm(form,data=dat,family=binomial(),control=glm.control(maxit=50L)),
    warning=function(w){warn <<- c(warn,conditionMessage(w));invokeRestart("muffleWarning")})
  cf <- coef(fit); p <- fitted(fit)
  if(length(warn)||!fit$converged||fit$rank!=length(cf)||any(!is.finite(cf))||
     any(!is.finite(p))||any(p<=1e-7|p>=1-1e-7))
    stop(paste(c("fit_gate_failure",warn),collapse=";"))
  fit
}
estimate <- function(dat,fit) {
  out <- numeric()
  for(term in c("vent","cvvh")) {
    zero <- one <- dat; zero[[term]]<-0L; one[[term]]<-1L
    r0 <- mean(predict(fit,newdata=zero,type="response"))
    r1 <- mean(predict(fit,newdata=one,type="response"))
    out[paste0(term,"_risk_no")]<-r0
    out[paste0(term,"_risk_yes")]<-r1
    out[paste0(term,"_rd_pp")]<-100*(r1-r0)
    out[paste0(term,"_or")]<-exp(coef(fit)[[term]])
  }
  stopifnot(all(is.finite(out))); out
}
fits <- lapply(forms,function(f)fit_checked(data,f))
points <- lapply(fits,function(f)estimate(data,f))
clusters <- split(seq_len(nrow(data)),data$patientid)
ids <- names(clusters); B <- 1000L; set.seed(20260930L)
boot <- array(NA_real_,dim=c(B,8L,2L),dimnames=list(NULL,names(points[[1]]),names(forms)))
coeff <- lapply(fits,function(f)matrix(NA_real_,nrow=B,ncol=length(coef(f)),dimnames=list(NULL,names(coef(f)))))
membership <- vector("list",B); failures <- list()
for(b in seq_len(B)) {
  drawn <- sample(ids,length(ids),replace=TRUE)
  index <- unlist(clusters[drawn],use.names=FALSE)
  membership[[b]] <- index
  sample_data <- data[index,,drop=FALSE]
  for(m in names(forms)) {
    tryCatch({f<-fit_checked(sample_data,forms[[m]])
      boot[b,,m]<-estimate(sample_data,f); coeff[[m]][b,]<-coef(f)},
      error=function(e){failures[[length(failures)+1L]] <<-
        data.frame(draw=b,model=m,reason=conditionMessage(e))})
  }
  if(b%%100L==0L)cat("paired draws completed: ",b,"/",B,"\n",sep="")
}
valid <- vapply(seq_len(B),function(b)all(is.finite(boot[b,,])),logical(1))
can_ci <- sum(valid)>=990L
ci <- function(x)if(can_ci)as.numeric(quantile(x[valid],c(.025,.975),type=7))else c(NA_real_,NA_real_)
rows <- list(); k<-0L
for(m in names(forms))for(metric in names(points[[m]])) {
  lim<-ci(boot[,metric,m]); k<-k+1L
  rows[[k]]<-data.frame(model=m,metric=metric,estimate=unname(points[[m]][metric]),
    ci_low=lim[1],ci_high=lim[2],n_events=634L,n_patients=624L,n_restarts=181L,
    paired_bootstrap_success=sum(valid),paired_bootstrap_failure=sum(!valid))
}
write.csv(do.call(rbind,rows),outputs[1],row.names=FALSE,na="")
change_rows<-lapply(c("vent","cvvh"),function(term){metric<-paste0(term,"_rd_pp")
  delta<-boot[,metric,"with_dose_band"]-boot[,metric,"without_dose_band"]; lim<-ci(delta)
  data.frame(indicator=term,contrast="with_minus_without_dose_band",
    rd_change_pp=unname(points$with_dose_band[metric]-points$without_dose_band[metric]),
    ci_low=lim[1],ci_high=lim[2],paired_bootstrap_success=sum(valid))})
write.csv(do.call(rbind,change_rows),outputs[2],row.names=FALSE,na="")
replicates<-data.frame(draw=seq_len(B),jointly_valid=valid)
for(m in names(forms))for(metric in names(points[[m]]))
  replicates[[paste(m,metric,sep="__")]]<-boot[,metric,m]
write.csv(replicates,outputs[3],row.names=FALSE,na="")
coef_rows<-list(); k<-0L
for(m in names(forms))for(term in names(coef(fits[[m]]))){lim<-ci(coeff[[m]][,term]);k<-k+1L
  coef_rows[[k]]<-data.frame(model=m,term=term,log_odds_coefficient=unname(coef(fits[[m]])[term]),
    ci_low=lim[1],ci_high=lim[2])}
write.csv(do.call(rbind,coef_rows),outputs[4],row.names=FALSE,na="")
fail_table<-if(length(failures))do.call(rbind,failures)else
  data.frame(draw=integer(),model=character(),reason=character())
write.csv(fail_table,outputs[5],row.names=FALSE)
saveRDS(list(data=data,membership=membership,forms=forms,seed=20260930L,points=points,
  boot=boot,coefficient_boot=coeff),file.path(stage,"same634_reproduction.rds"),compress="gzip")
diagnostics<-lapply(fits,function(f){inf<-cooks.distance(f)
  list(formula=paste(deparse(formula(f)),collapse=" "),converged=f$converged,rank=f$rank,
    parameter_count=length(coef(f)),fitted_probability_range=range(fitted(f)),
    largest_10_cook_distances=sort(inf,decreasing=TRUE)[seq_len(10)],
    largest_10_abs_dfbeta=sort(apply(abs(dfbetas(f)),1,max),decreasing=TRUE)[seq_len(10)])})
audit<-list(started=as.character(started),finished=as.character(Sys.time()),
  elapsed_seconds=as.numeric(difftime(Sys.time(),started,units="secs")),R_version=R.version.string,
  seed=20260930L,bootstrap_draws=1000L,patient_clusters=624L,events=634L,restarts=181L,
  jointly_successful=sum(valid),jointly_failed=sum(!valid),same_rows_both_models=TRUE,
  standardisation="each bootstrap draw's sampled events",percentile_type=7,
  age_counts=as.list(table(data$age)),sex_counts=as.list(table(data$sex)),
  peak_band_counts=as.list(table(data$peak_band)),joint_support_cells=unclass(table(data$vent,data$cvvh)),
  diagnostics=diagnostics,no_individual_predictions_exported=TRUE)
jsonlite::write_json(audit,outputs[6],auto_unbox=TRUE,pretty=TRUE,digits=NA)
cat("Completed same-634 paired comparison. Successful paired draws: ",sum(valid),"\n",sep="")
