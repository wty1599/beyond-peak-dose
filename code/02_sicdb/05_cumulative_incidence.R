options(stringsAsFactors=FALSE,warn=1,digits=17)
suppressPackageStartupMessages({library(data.table);library(jsonlite);library(digest);library(survival)})
setDTthreads(1L);stage<-commandArgs(TRUE)[1];stopifnot(stage%in%c('S2','S3','S4'))
source('code/00_setup/paths.R');out<-private_output('sicdb_clinical');start<-Sys.time()
src<-release_path('sicdb_clinical_input');before<-digest(file=src,algo='sha256');d<-fread(src)
source('code/00_setup/aalen_johansen.R')
stopifnot(nrow(d)==2044,uniqueN(d$PatientID)==1996,sum(d$restart_binary,na.rm=TRUE)==684)
d[,`:=`(MV=predictor_mechanical_ventilation_t0,CRRT=recorded_crrt_active_pre6h)]
if(stage=='S3')d<-d[!is.na(restart_binary)]
groups<-list();spec<-data.table()
add<-function(name,ix,landmark=4,times=c(12,24,72)){
 groups[[length(groups)+1L]]<<-ix
 spec<<-rbind(spec,data.table(group=name,landmark_hours=landmark,report_time_hours=list(times)))
}
if(stage=='S2'){
 add('Overall',seq_len(nrow(d)),times=c(12,24,48,72))
 for(i in 1:3){v<-d$predictor_episode_peak_nee;ix<-which(if(i==1)v<.2 else if(i==2)v>=.2&v<.5 else v>=.5);add(c('<0.20','0.20-<0.50','>=0.50')[i],ix,times=c(12,24,48,72))}
}
states<-c('RESTART','DEATH','NORMAL_CARE_END','ADMIN_HORIZON','CARE_END_UNCLASSIFIED')
describe<-function(z,group,joint=FALSE){data.table(group=group,cases=nrow(z),patients=uniqueN(z$PatientID),restart_n=sum(z$restart_binary==1,na.rm=TRUE),restart_known=sum(!is.na(z$restart_binary)),restart_unknown=sum(is.na(z$restart_binary)),restart_proportion=if(!joint)mean(z$restart_binary,na.rm=TRUE)else NA_real_,death_n=sum(z$death_72h_binary==1,na.rm=TRUE),death_known=sum(!is.na(z$death_72h_binary)),death_unknown=sum(is.na(z$death_72h_binary)),death_proportion=if(!joint)mean(z$death_72h_binary,na.rm=TRUE)else NA_real_,first_restart=sum(z$cif_state=='RESTART'),first_death=sum(z$cif_state=='DEATH'),first_normal_case_end=sum(z$cif_state=='NORMAL_CARE_END'),first_horizon=sum(z$cif_state=='ADMIN_HORIZON'),first_unclassified_end=sum(z$cif_state=='CARE_END_UNCLASSIFIED'),joint_counts_only=joint)}
if(stage=='S3'){
 av<-rbindlist(lapply(c('MV','CRRT'),function(k)data.table(variable=k,eligible_cases=nrow(d),known=sum(!is.na(d[[k]])),unknown=sum(is.na(d[[k]])),recorded_or_proxy_positive=sum(d[[k]]==1,na.rm=TRUE),recorded_or_proxy_zero=sum(d[[k]]==0,na.rm=TRUE))))
 fwrite(av,file.path(out,'S3_support_availability.csv'));obs<-list()
 for(k in c('MV','CRRT'))for(v in 0:1){ix<-which(d[[k]]==v);nm<-paste0(k,'=',v);add(nm,ix);obs[[length(obs)+1L]]<-describe(d[ix],nm)}
 for(m in 0:1)for(r in 0:1)obs[[length(obs)+1L]]<-describe(d[MV==m & CRRT==r],paste0('MV=',m,';CRRT=',r),TRUE)
 desc<-rbindlist(obs);stopifnot(all(desc$cases==rowSums(desc[,.(first_restart,first_death,first_normal_case_end,first_horizon,first_unclassified_end)])))
 fwrite(desc,file.path(out,'S3_organ_support_stratified.csv'))
}
if(stage=='S4'){
 flow<-list()
 for(lm in c(12,24))for(g in c('Overall','MV=0','MV=1')){
  ix<-if(g=='Overall')which(d$first_time_hours>lm)else which(d$first_time_hours>lm & !is.na(d$MV)&d$MV==as.integer(sub('MV=','',g)))
  add(g,ix,lm,72);z<-d[ix];rr<-describe(z,g);rr[,`:=`(landmark_hours=lm,analysis_end_hours=72,MV_unknown=sum(is.na(z$MV)))];flow[[length(flow)+1L]]<-rr
 }
 fwrite(rbindlist(flow),file.path(out,'S4_conditional_population_and_first_events.csv'))
}
fail<-data.table(analysis=character(),bootstrap=integer(),group=character(),reason=character());fwrite(fail,file.path(out,paste0(stage,'_bootstrap_failures.csv')))
uid<-unique(d$PatientID);pid<-match(d$PatientID,uid);points<-curves<-boots<-checks<-list()
for(an in c('NORMAL_CASE_END_COMPETING','NORMAL_CASE_END_CENSORED')){
 status<-if(an=='NORMAL_CASE_END_COMPETING')d$first_status_case_competing else d$first_status_main
 pp<-list();positions<-list();k<-0L
 for(gi in seq_along(groups)){
  ix<-groups[[gi]];tt<-spec$report_time_hours[[gi]];lm<-spec$landmark_hours[gi];stopifnot(length(ix)>0)
  z<-aj(d$first_time_hours[ix],status[ix]);vv<-stepat(z,tt,'restart');ref<-aj_reference(d$first_time_hours[ix],status[ix],tt)
  sf<-survfit(Surv(d$first_time_hours[ix]-lm,factor(status[ix],levels=0:3,labels=c('censor','restart','death','care_end')))~1)
  delta<-max(abs(sf$pstate[,match('restart',sf$states)]-stepat(z,sf$time+lm,'restart')))
  stopifnot(max(abs(vv-ref))<1e-12,delta<1e-10)
  checks[[length(checks)+1L]]<-data.table(analysis=an,group=spec$group[gi],landmark_hours=lm,loop_delta=max(abs(vv-ref)),survfit_delta=delta)
  nrs<-vapply(tt,function(t)sum(d$first_time_hours[ix]>=t),integer(1))
  pp[[gi]]<-data.table(analysis=an,group=spec$group[gi],landmark_hours=lm,time_hours=tt,cases=length(ix),patients=uniqueN(d$PatientID[ix]),subsequent_restarts=sum(status[ix]==1),estimate=vv,at_risk_before=nrs)
  positions[[gi]]<-k+seq_along(tt);k<-k+length(tt)
  z<-rbind(data.table(time_hours=lm,restart=0,death=0,care_end=0,event_free=1,at_risk_before=length(ix)),z)
  z[,`:=`(analysis=an,group=spec$group[gi],landmark_hours=lm)];curves[[length(curves)+1L]]<-z
 }
 point<-rbindlist(pp);bt<-riskbt<-matrix(NA_real_,1000,nrow(point));set.seed(202609094L)
 for(b in 1:1000){
  wt<-tabulate(sample.int(length(uid),length(uid),replace=TRUE),nbins=length(uid))[pid]
  for(gi in seq_along(groups)){
   ix<-groups[[gi]];tt<-spec$report_time_hours[[gi]]
   result<-tryCatch({stopifnot(sum(wt[ix])>0);z<-aj(d$first_time_hours[ix],status[ix],wt[ix]);v<-stepat(z,tt,'restart');stopifnot(all(is.finite(v)));list(v=v,n=vapply(tt,function(t)sum(wt[ix][d$first_time_hours[ix]>=t]),numeric(1)))},error=function(e)e)
   if(inherits(result,'error')){fail<-rbind(fail,data.table(analysis=an,bootstrap=b,group=spec$group[gi],reason=conditionMessage(result)));fwrite(fail,file.path(out,paste0(stage,'_bootstrap_failures.csv')));stop('CIF bootstrap failure; retained reason; no fallback')}
   bt[b,positions[[gi]]]<-result$v;riskbt[b,positions[[gi]]]<-result$n
  }
  if(b%%250==0){cat(stage,an,b,'/1000\n');flush.console()}
 }
 point[,`:=`(lower95=apply(bt,2,quantile,.025,type=7),upper95=apply(bt,2,quantile,.975,type=7),bootstrap_success=1000L,seed=202609094L,min_bootstrap_at_risk=apply(riskbt,2,min),bootstrap_zero_risk_set=colSums(riskbt==0),source_operation='new calculation')]
 stopifnot(all(is.finite(point$lower95)),all(is.finite(point$upper95)))
 if(any(point$at_risk_before==0)){fwrite(point,file.path(out,paste0(stage,'_risk_set_stop.csv')));stop('Observed risk set exhausted at reporting time; author review required')}
 points[[length(points)+1L]]<-point
 boots[[length(boots)+1L]]<-data.table(analysis=an,bootstrap=rep(1:1000,nrow(point)),group=rep(point$group,each=1000),landmark_hours=rep(point$landmark_hours,each=1000),time_hours=rep(point$time_hours,each=1000),estimate=as.vector(bt),at_risk_before=as.vector(riskbt))
}
res<-rbindlist(points)
stem<-c(S2='S2_cif_overall_and_peak_strata',S3='S3_organ_support_cif',S4='S4_conditional_residual_risk')[stage]
fwrite(res,file.path(out,paste0(stem,'.csv')));fwrite(rbindlist(curves),file.path(out,paste0(stage,'_full_curves.csv')));fwrite(rbindlist(boots),file.path(out,paste0(stage,'_patient_cluster_bootstrap.csv')));fwrite(rbindlist(checks),file.path(out,paste0(stage,'_AJ_crosschecks.csv')))
stopifnot(digest(file=src,algo='sha256')==before)
write_json(list(status='PASS_PENDING_INDEPENDENT_REVIEW',base_cases=nrow(d),base_patients=uniqueN(d$PatientID),point_rows=nrow(res),bootstrap_replicates_each_version=1000,seed=202609094,bootstrap_failures=nrow(fail),minimum_risk_set=min(res$at_risk_before),minimum_bootstrap_risk_set=min(res$min_bootstrap_at_risk),source_unchanged=TRUE,patient_cluster_sampling='Full stage mother-cohort patient multiplicity carried to prespecified strata and conditional risk sets',elapsed_seconds=as.numeric(difftime(Sys.time(),start,units='secs'))),file.path(out,paste0(stage,'_QC.json')),auto_unbox=TRUE,pretty=TRUE)
print(res)
