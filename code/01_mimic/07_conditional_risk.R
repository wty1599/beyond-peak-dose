options(stringsAsFactors=FALSE,warn=1,digits=17)
Sys.setenv(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
suppressPackageStartupMessages({library(data.table);library(survival);library(jsonlite);library(digest)})
setDTthreads(1L)
source('code/00_setup/paths.R')
root<-private_output('mimic_conditional')
started<-Sys.time()
paths<-c(cif_input=release_path('mimic_clinical_input'),predictors=release_path('mimic_observed_input'),helper='code/00_setup/aalen_johansen.R')
hashes<-vapply(paths,function(p)digest(file=p,algo='sha256'),character(1))
source(paths['helper'])
d<-fread(paths['cif_input'])
p<-fread(paths['predictors'],select=c('stay_id','subject_id','predictor_mechanical_ventilation_t0'))
stopifnot(nrow(d)==1958L,uniqueN(d$subject_id)==1856L,uniqueN(d$stay_id)==1958L,
 nrow(p)==1958L,uniqueN(p$stay_id)==1958L,all(d$entry_hours==4),
 all(d$time_hours>=4 & d$time_hours<=72),all(d$cohort=='MIMIC'),all(d$version=='minimal_1958'),
 sum(d$state=='RESTART')==680L,all(d$state%in%c('RESTART','DEATH','NORMAL_CARE_END','ADMIN_HORIZON')),
 all(d$status_s20==fcase(d$state=='RESTART',1L,d$state=='DEATH',2L,d$state=='NORMAL_CARE_END',3L,default=0L)),
 all(d$status_main==fcase(d$state=='RESTART',1L,d$state=='DEATH',2L,default=0L)))
j<-match(d$stay_id,p$stay_id)
stopifnot(!anyNA(j),all(d$subject_id==p$subject_id[j]))
d[,MV:=p$predictor_mechanical_ventilation_t0[j]]
stopifnot(!anyNA(d$MV),all(d$MV%in%0:1))
spec<-CJ(landmark_hours=c(12,24),group=c('Overall','MV=0','MV=1'),sorted=FALSE)
getbase<-function(g)if(g=='Overall')seq_len(nrow(d))else which(d$MV==as.integer(sub('MV=','',g)))
indices<-lapply(seq_len(nrow(spec)),function(i){ix<-getbase(spec$group[i]);ix[d$time_hours[ix]>spec$landmark_hours[i]]})
flow<-rbindlist(lapply(seq_len(nrow(spec)),function(i){
 base<-getbase(spec$group[i]);ix<-indices[[i]];s<-spec$landmark_hours[i];z<-d[ix];a<-d[base]
 data.table(landmark_hours=s,remaining_hours=72-s,group=spec$group[i],baseline_events=nrow(a),
 eligible_events=nrow(z),eligible_patients=uniqueN(z$subject_id),
 excluded_restart_by_landmark=sum(a$state=='RESTART' & a$time_hours<=s),
 excluded_death_by_landmark=sum(a$state=='DEATH' & a$time_hours<=s),
 excluded_normal_end_by_landmark=sum(a$state=='NORMAL_CARE_END' & a$time_hours<=s),
 exact_landmark_restart=sum(a$state=='RESTART' & a$time_hours==s),
 exact_landmark_death=sum(a$state=='DEATH' & a$time_hours==s),
 exact_landmark_normal_end=sum(a$state=='NORMAL_CARE_END' & a$time_hours==s),
 remaining_first_restart=sum(z$state=='RESTART'),remaining_first_death=sum(z$state=='DEATH'),
 remaining_first_normal_end=sum(z$state=='NORMAL_CARE_END'),remaining_admin_horizon=sum(z$state=='ADMIN_HORIZON'))
}))
stopifnot(all(flow$eligible_events>0),all(flow$baseline_events==flow$eligible_events+flow$excluded_restart_by_landmark+flow$excluded_death_by_landmark+flow$excluded_normal_end_by_landmark),
 all(flow$eligible_events==flow$remaining_first_restart+flow$remaining_first_death+flow$remaining_first_normal_end+flow$remaining_admin_horizon))
fwrite(flow,file.path(root,'conditional_population_and_first_events.csv'))
fail<-data.table(analysis=character(),bootstrap=integer(),landmark_hours=numeric(),group=character(),reason=character())
fwrite(fail,file.path(root,'bootstrap_failures.csv'))
uid<-unique(d$subject_id);pid<-match(d$subject_id,uid)
allpoints<-allboots<-checks<-list()
for(an in c('S20_NORMAL_END_COMPETING','PRIMARY_CARE_END_CENSORED')){
 status<-if(an=='S20_NORMAL_END_COMPETING')d$status_s20 else d$status_main
 seed<-if(an=='S20_NORMAL_END_COMPETING')2026072941L else 2026072906L
 point<-copy(flow)
 estimates<-vapply(seq_len(nrow(spec)),function(i){
  ix<-indices[[i]];s<-spec$landmark_hours[i];z<-aj(d$time_hours[ix]-s,status[ix]);est<-stepat(z,72-s,'restart')
  ref<-aj_reference(d$time_hours[ix]-s,status[ix],72-s)
  sf<-survfit(Surv(d$time_hours[ix]-s,factor(status[ix],levels=0:3,labels=c('censor','restart','death','care_end')))~1)
  sfval<-tail(sf$pstate[,match('restart',sf$states)],1)
  checks[[length(checks)+1L]]<<-data.table(analysis=an,landmark_hours=s,group=spec$group[i],
    event_loop_difference=abs(est-ref),survfit_difference=abs(est-sfval),
    count_ratio_difference=if(an=='S20_NORMAL_END_COMPETING')abs(est-mean(d$state[ix]=='RESTART'))else NA_real_)
  stopifnot(abs(est-ref)<1e-12,abs(est-sfval)<1e-10)
  if(an=='S20_NORMAL_END_COMPETING')stopifnot(abs(est-mean(d$state[ix]=='RESTART'))<1e-12)
  est
 },numeric(1))
 set.seed(seed)
 bt<-matrix(NA_real_,1000L,nrow(spec));nr<-matrix(NA_real_,1000L,nrow(spec))
 for(b in 1:1000){
  weights<-tabulate(sample.int(length(uid),length(uid),replace=TRUE),nbins=length(uid))[pid]
  for(i in seq_len(nrow(spec))){
   ix<-indices[[i]];s<-spec$landmark_hours[i]
   result<-tryCatch({
    nn<-sum(weights[ix]);if(nn<=0)stop('Empty conditional stratum')
    z<-aj(d$time_hours[ix]-s,status[ix],weights[ix]);v<-stepat(z,72-s,'restart')
    if(!is.finite(v)||v<0||v>1)stop('Invalid conditional CIF')
    list(value=v,events=nn)
   },error=function(e)e)
   if(inherits(result,'error')){
    fail<-rbind(fail,data.table(analysis=an,bootstrap=b,landmark_hours=s,group=spec$group[i],reason=conditionMessage(result)))
    fwrite(fail,file.path(root,'bootstrap_failures.csv'))
   }else{bt[b,i]<-result$value;nr[b,i]<-result$events}
  }
  if(b%%250L==0){cat(an,' ',b,'/1000\n',sep='');flush.console()}
 }
 point[,`:=`(analysis=an,estimate=estimates,lower95=vapply(seq_len(ncol(bt)),function(j)if(all(is.finite(bt[,j])))quantile(bt[,j],.025,type=7)else NA_real_,numeric(1)),
  upper95=vapply(seq_len(ncol(bt)),function(j)if(all(is.finite(bt[,j])))quantile(bt[,j],.975,type=7)else NA_real_,numeric(1)),
  bootstrap_success=colSums(is.finite(bt)),bootstrap_target=1000L,seed=seed,
  min_bootstrap_eligible_events=apply(nr,2,min,na.rm=TRUE),
  presentation=if(an=='S20_NORMAL_END_COMPETING')'v3.1_main_presentation'else'v3.1_supplement_original_SAP_primary')]
 point[,interval_status:=ifelse(bootstrap_success==1000,'ESTIMATED','NOT_ESTIMABLE_FAILED_REPLICATES')]
 allpoints[[length(allpoints)+1L]]<-point
 allboots[[length(allboots)+1L]]<-data.table(analysis=an,bootstrap=rep(1:1000,ncol(bt)),
  landmark_hours=rep(spec$landmark_hours,each=1000),group=rep(spec$group,each=1000),
  estimate=as.vector(bt),sampled_eligible_events=as.vector(nr))
}
out<-rbindlist(allpoints);btout<-rbindlist(allboots);check<-rbindlist(checks)
fwrite(out,file.path(root,'conditional_restart_risk.csv'))
fwrite(btout,file.path(root,'conditional_bootstrap_aggregate.csv'))
fwrite(check,file.path(root,'AJ_method_crosschecks.csv'))
after<-vapply(paths,function(p)digest(file=p,algo='sha256'),character(1))
fwrite(data.table(role=names(paths),path=unname(paths),sha256_before=unname(hashes),sha256_after=unname(after),unchanged=unname(hashes==after)),file.path(root,'source_hash_comparison.csv'))
stopifnot(all(after==hashes),nrow(fail)==0L,all(out$bootstrap_success==1000L),nrow(out)==12L)
csvs<-list.files(root,pattern='\\.csv$',full.names=TRUE)
identifier_check<-rbindlist(lapply(csvs,function(f){cols<-names(fread(f,nrows=0));bad<-intersect(tolower(cols),c('subject_id','hadm_id','stay_id','case_id','admission_id','patient_id','prediction','predicted_risk'))
 data.table(file=basename(f),forbidden_identifier_columns=length(bad))}))
fwrite(identifier_check,file.path(root,'aggregate_identifier_check.csv'))
stopifnot(all(identifier_check$forbidden_identifier_columns==0))
qc<-list(PASS=TRUE,cohort_events=1958L,cohort_patients=1856L,restarts=680L,
 landmark_hours=c(12L,24L),horizon_from_t0_hours=72L,eligible_definition='first_event_time strictly greater than landmark',
 groups=c('Overall','MV=0','MV=1'),MV_measured_at='t0',bootstrap_success_per_cell=1000L,result_cells=12L,
 failures=0L,seed_S20=2026072941L,seed_original_primary=2026072906L,
 sampling='Original full-cohort patient multiplicities; all events carried jointly; landmark filtering within each draw',
 independent_QC_status='PENDING',global_protection_after='PENDING',source_hashes_unchanged=TRUE,
 model_fits=0L,event_readjudications=0L,external_analyses=0L,group_tests=0L,row_level_outputs=0L,
 elapsed_seconds=as.numeric(difftime(Sys.time(),started,units='secs')))
write_json(qc,file.path(root,'CONDITIONAL_QC.json'),pretty=TRUE,auto_unbox=TRUE,digits=17)
capture.output(sessionInfo(),file=file.path(root,'R_session_info.txt'))
print(out[,.(analysis,landmark_hours,group,eligible_events,eligible_patients,remaining_first_restart,estimate,lower95,upper95)])
print(qc)
