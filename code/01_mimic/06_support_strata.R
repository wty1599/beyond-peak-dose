options(stringsAsFactors=FALSE, warn=1, digits=17)
Sys.setenv(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
suppressPackageStartupMessages({library(data.table); library(survival); library(jsonlite); library(digest)})
setDTthreads(1L)
source('code/00_setup/paths.R')
root<-private_output('mimic_support')
started<-Sys.time()
paths<-c(cohort=release_path('mimic_selected_cohort'),predictors=release_path('mimic_observed_input'),cif_input=release_path('mimic_clinical_input'),aj_helper='code/00_setup/aalen_johansen.R')
input_hash<-vapply(paths,function(p)digest(file=p,algo='sha256'),character(1))
source(paths['aj_helper'])
d <- fread(paths['cif_input'])
s <- fread(paths['cohort'],select=c('stay_id','subject_id','restart_72h','death_72h','cif_status'))
p <- fread(paths['predictors'],select=c('stay_id','subject_id','predictor_mechanical_ventilation_t0','predictor_rrt_t0'))
stopifnot(nrow(d)==1958L, uniqueN(d$subject_id)==1856L, uniqueN(d$stay_id)==1958L,
          nrow(s)==1958L,nrow(p)==1958L, uniqueN(s$stay_id)==1958L,uniqueN(p$stay_id)==1958L)
si <- match(d$stay_id,s$stay_id); pi <- match(d$stay_id,p$stay_id)
stopifnot(!anyNA(si),!anyNA(pi),all(d$subject_id==s$subject_id[si]),all(d$subject_id==p$subject_id[pi]),
          all(d$state==s$cif_status[si]),all(d$time_hours>=4 & d$time_hours<=72),all(d$entry_hours==4))
d[,`:=`(MV=p$predictor_mechanical_ventilation_t0[pi], RRT=p$predictor_rrt_t0[pi],
        restart_72h=s$restart_72h[si],death_72h=s$death_72h[si])]
stopifnot(!anyNA(d$MV),!anyNA(d$RRT),all(d$MV%in%0:1),all(d$RRT%in%0:1),
          !anyNA(d$restart_72h),sum(d$restart_72h)==680L,
          sum(d$death_72h==1,na.rm=TRUE)==74L,sum(is.na(d$death_72h))==16L,
          all(d$status_main==fifelse(d$state=='RESTART',1L,fifelse(d$state=='DEATH',2L,0L))),
          all(d$status_s20==fifelse(d$state=='RESTART',1L,fifelse(d$state=='DEATH',2L,fifelse(d$state=='NORMAL_CARE_END',3L,0L)))))
summary_rows <- function(ix, group_type, group_label, mv=NA_integer_,rrt=NA_integer_) {
 z<-d[ix]
 data.table(group_type=group_type, group=group_label,MV_value=mv,RRT_value=rrt,
  events=nrow(z),patients=uniqueN(z$subject_id),
  restart_n=sum(z$restart_72h==1,na.rm=TRUE),restart_evaluable_n=sum(!is.na(z$restart_72h)),
  restart_unknown_n=sum(is.na(z$restart_72h)),restart_proportion=mean(z$restart_72h,na.rm=TRUE),
  death_n=sum(z$death_72h==1,na.rm=TRUE),death_evaluable_n=sum(!is.na(z$death_72h)),
  death_unknown_n=sum(is.na(z$death_72h)),death_proportion=mean(z$death_72h,na.rm=TRUE),
  first_restart=sum(z$state=='RESTART'),first_death=sum(z$state=='DEATH'),
  first_normal_care_end=sum(z$state=='NORMAL_CARE_END'),first_72h_horizon=sum(z$state=='ADMIN_HORIZON'))
}
groups <- data.table(variable=rep(c('MV','RRT'),each=2),value=rep(0:1,2))
descriptive <- rbindlist(c(list(summary_rows(seq_len(nrow(d)),'Overall','Overall')),
 lapply(seq_len(nrow(groups)),function(i) { v<-groups$variable[i];a<-groups$value[i];
 summary_rows(which(d[[v]]==a),v,paste0(v,'=',a),if(v=='MV')a else NA_integer_,if(v=='RRT')a else NA_integer_) }),
 lapply(0:3,function(i){mv<-i%/%2;rrt<-i%%2;summary_rows(which(d$MV==mv & d$RRT==rrt),'MV_x_RRT',paste0('MV=',mv,';RRT=',rrt),mv,rrt)})))
stopifnot(all(descriptive$events==descriptive$first_restart+descriptive$first_death+descriptive$first_normal_care_end+descriptive$first_72h_horizon))
descriptive[group_type=='MV_x_RRT',`:=`(restart_proportion=NA_real_,death_proportion=NA_real_)]
fwrite(descriptive,file.path(root,'organ_support_observed_outcomes.csv'))

times <- c(12,24,72)
uid<-unique(d$subject_id); pid<-match(d$subject_id,uid)
failures<-data.table(analysis=character(),bootstrap=integer(),variable=character(),value=integer(),reason=character())
fwrite(failures,file.path(root,'bootstrap_failures.csv'))
allpoints<-allboots<-riskrows<-checks<-list()
for(an in c('PRIMARY_CARE_END_CENSORED','S20_NORMAL_END_COMPETING')){
  status<-if(an=='PRIMARY_CARE_END_CENSORED')d$status_main else d$status_s20
  seed<-if(an=='PRIMARY_CARE_END_CENSORED')2026072906L else 2026072941L
  point<-rbindlist(lapply(seq_len(nrow(groups)),function(gi){
    v<-groups$variable[gi];val<-groups$value[gi];ix<-which(d[[v]]==val)
    z<-aj(d$time_hours[ix],status[ix]);est<-stepat(z,times,'restart')
    ref<-aj_reference(d$time_hours[ix],status[ix],times,cause=1L)
    sf<-survfit(Surv(d$time_hours[ix]-4,factor(status[ix],levels=0:3,labels=c('censor','restart','death','care_end')))~1)
    sfdelta<-max(abs(sf$pstate[,match('restart',sf$states)]-stepat(z,sf$time+4,'restart')))
    checks[[length(checks)+1L]]<<-data.table(analysis=an,variable=v,value=val,reference='independent_event_loop',max_absolute_difference=max(abs(est-ref)))
    checks[[length(checks)+1L]]<<-data.table(analysis=an,variable=v,value=val,reference='survival_survfit',max_absolute_difference=sfdelta)
    stopifnot(max(abs(est-ref))<1e-12,sfdelta<1e-10)
    nr<-vapply(times,function(t)sum(d$time_hours[ix]>=t),integer(1))
    riskrows[[length(riskrows)+1L]]<<-data.table(analysis=an,variable=v,value=val,time_hours=times,at_risk_before=nr)
    data.table(analysis=an,variable=v,value=val,time_hours=times,estimate=est,events=length(ix),patients=uniqueN(d$subject_id[ix]),at_risk_before=nr)
  }))
  set.seed(seed)
  bt<-matrix(NA_real_,nrow=1000L,ncol=nrow(point))
  riskbt<-matrix(NA_real_,nrow=1000L,ncol=nrow(point))
  for(b in 1:1000){
    frequencies<-tabulate(sample.int(length(uid),length(uid),replace=TRUE),nbins=length(uid))
    weights<-frequencies[pid]
    for(gi in seq_len(nrow(groups))){
      v<-groups$variable[gi];val<-groups$value[gi];ix<-which(d[[v]]==val);cols<-(gi-1L)*3L+1:3
      result<-tryCatch({
        if(sum(weights[ix])<=0)stop('No sampled events in the prespecified stratum')
        zz<-aj(d$time_hours[ix],status[ix],weights[ix]);out<-stepat(zz,times,'restart')
        if(any(!is.finite(out)))stop('Nonfinite CIF')
        list(out=out,nr=vapply(times,function(t)sum(weights[ix][d$time_hours[ix]>=t]),numeric(1)))
      },error=function(e)e)
      if(inherits(result,'error')){
        failures<-rbind(failures,data.table(analysis=an,bootstrap=b,variable=v,value=val,reason=conditionMessage(result)))
        fwrite(failures,file.path(root,'bootstrap_failures.csv'))
      }else{bt[b,cols]<-result$out;riskbt[b,cols]<-result$nr}
    }
    if(b%%250L==0){cat(an,' bootstrap ',b,'/1000\n',sep='');flush.console()}
  }
  point[,`:=`(lower95=vapply(seq_len(ncol(bt)),function(j)if(all(is.finite(bt[,j])))quantile(bt[,j],.025,type=7)else NA_real_,numeric(1)),
              upper95=vapply(seq_len(ncol(bt)),function(j)if(all(is.finite(bt[,j])))quantile(bt[,j],.975,type=7)else NA_real_,numeric(1)),
              bootstrap_success=colSums(is.finite(bt)),bootstrap_target=1000L,seed=seed,
              min_bootstrap_at_risk=apply(riskbt,2,min,na.rm=TRUE),
              bootstrap_zero_risk_set=colSums(riskbt==0,na.rm=TRUE))]
  point[,interval_status:=ifelse(bootstrap_success==1000L,'ESTIMATED','NOT_ESTIMABLE_FAILED_REPLICATES')]
  allpoints[[length(allpoints)+1L]]<-point
  allboots[[length(allboots)+1L]]<-data.table(analysis=an,bootstrap=rep(1:1000,ncol(bt)),
    variable=rep(point$variable,each=1000),value=rep(point$value,each=1000),
    time_hours=rep(point$time_hours,each=1000),estimate=as.vector(bt),at_risk_before=as.vector(riskbt))
}
points<-rbindlist(allpoints);boots<-rbindlist(allboots);check<-rbindlist(checks)
fwrite(points,file.path(root,'organ_support_restart_cif_12_24_72h.csv'))
fwrite(boots,file.path(root,'patient_cluster_bootstrap_aggregate.csv'))
fwrite(rbindlist(riskrows),file.path(root,'risk_sets.csv'))
fwrite(check,file.path(root,'AJ_method_crosschecks.csv'))
after<-vapply(paths,function(p)digest(file=p,algo='sha256'),character(1))
unchanged<-after==input_hash
fwrite(data.table(role=names(paths),path=unname(paths),sha256_before=unname(input_hash),sha256_after=unname(after),unchanged=unname(unchanged)),file.path(root,'source_hash_comparison.csv'))
stopifnot(all(unchanged),nrow(failures)==0L,all(points$bootstrap_success==1000L),all(points$at_risk_before>0),all(points$bootstrap_zero_risk_set==0L))
qc<-list(PASS=TRUE,events=nrow(d),patients=uniqueN(d$subject_id),restarts=sum(d$restart_72h),death_evaluable=sum(!is.na(d$death_72h)),
 deaths=sum(d$death_72h,na.rm=TRUE),death_unknown=sum(is.na(d$death_72h)),
 primary_seed=2026072906L,S20_seed=2026072941L,bootstrap_draws_per_analysis=1000L,patient_cluster_sampling='global cohort patient multiplicity; all events carried together',
 stratification_variables=c('MV','RRT'),joint_states_counts_only=TRUE,cif_estimates=nrow(points),bootstrap_failures=nrow(failures),
 minimum_observed_at_risk=min(points$at_risk_before),minimum_bootstrap_at_risk=min(points$min_bootstrap_at_risk),
 unestimable_intervals=sum(points$interval_status!='ESTIMATED'),source_hashes_unchanged=all(unchanged),
 event_readjudications=0L,regression_models=0L,group_tests=0L,interaction_models=0L,
 point_reference_checks=nrow(check),reference_max_difference=max(check$max_absolute_difference),
 data_table_threads=getDTthreads(),elapsed_seconds=as.numeric(difftime(Sys.time(),started,units='secs')),
 independent_external_process_review='PENDING_ROOT_REVIEW',global_protected_after_check='PENDING_ROOT_FINALIZATION')
write_json(qc,file.path(root,'T1_QC.json'),pretty=TRUE,auto_unbox=TRUE,digits=17)
capture.output(sessionInfo(),file=file.path(root,'R_session_info.txt'))
print(descriptive);print(points);print(qc)
