options(stringsAsFactors=FALSE,warn=1)
Sys.setenv(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
suppressPackageStartupMessages({library(data.table);library(survival);library(jsonlite);library(digest)})
setDTthreads(1L)
source('code/00_setup/paths.R')
root<-private_output('mimic_clinical')
a<-file.path(root,'aggregate');w<-file.path(root,'restricted_work')
if(file.exists(file.path(a,'cif_qc.json')))stop('Completed CIF analysis is protected')
source('code/00_setup/aalen_johansen.R')
started<-Sys.time()
manifest<-fread(file.path(a,'source_manifest.csv'))
verify_sources<-function()stopifnot(all(vapply(manifest$path,function(p)digest(file=p,algo='sha256'),character(1))==manifest$sha256))
verify_sources()
input<-release_path('mimic_clinical_input');input_hash<-digest(file=input,algo='sha256')
d<-fread(input)
stopifnot(nrow(d)==1958L,uniqueN(d$subject_id)==1856L,uniqueN(d$stay_id)==1958L,
          all(d$entry_hours==4),all(d$time_hours>=4&d$time_hours<=72),all(d$cohort=='MIMIC'),all(d$version=='minimal_1958'))

tests<-list()
groups<-c('Overall','<0.2','0.2-<0.5','>=0.5');horizons<-c(12,24,72);causes<-c('restart','death','care_end')
grid<-seq(4,72,by=.25);curves<-estimates<-checks<-risks<-boots<-list();fullband<-NULL;fullboot<-NULL;band_delta<-0
uid<-unique(d$subject_id);pid<-match(d$subject_id,uid)
for(an in c('PRIMARY_CARE_END_CENSORED','S20_NORMAL_END_COMPETING')){
 st<-if(an=='PRIMARY_CARE_END_CENSORED')d$status_main else d$status_s20
 point<-list()
 for(g in groups){ix<-if(g=='Overall')seq_len(nrow(d)) else which(d$dose_group==g)
  z<-aj(d$time_hours[ix],st[ix]);z[,`:=`(cohort='MIMIC',version='minimal_1958',analysis=an,dose_group=g)]
  shifted<-d$time_hours[ix]-4
  sf<-survfit(Surv(shifted,factor(st[ix],levels=0:3,labels=c('censor','restart','death','care_end')))~1)
  for(cause in causes){
   delta<-max(abs(stepat(z,horizons,cause)-aj_reference(d$time_hours[ix],st[ix],horizons,cause=match(cause,causes))))
   checks[[length(checks)+1L]]<-data.table(analysis=an,dose_group=g,cause=cause,reference='independent_MIMIC_event_loop',max_absolute_difference=delta)
   stopifnot(delta<1e-12)
   delta_sf<-max(abs(sf$pstate[,match(cause,sf$states)]-stepat(z,sf$time+4,cause)))
   checks[[length(checks)+1L]]<-data.table(analysis=an,dose_group=g,cause=cause,reference='survival_survfit_all_observed_times',max_absolute_difference=delta_sf)
   stopifnot(delta_sf<1e-10)
   point[[length(point)+1L]]<-data.table(cohort='MIMIC',version='minimal_1958',analysis=an,dose_group=g,cause=cause,time_hours=horizons,
     estimate=stepat(z,horizons,cause),n=length(ix),patients=uniqueN(d$subject_id[ix]))
  }
  stopifnot(all(z$time_hours>=4),all(z$restart>=0),all(diff(z$restart)>=-1e-12),all(diff(z$death)>=-1e-12),all(diff(z$care_end)>=-1e-12),
            all(abs(z$restart+z$death+z$care_end+z$event_free-1)<1e-10))
  curves[[length(curves)+1L]]<-z
  rt<-c(4,12,24,48,72)
  risks[[length(risks)+1L]]<-data.table(cohort='MIMIC',version='minimal_1958',analysis=an,dose_group=g,time_hours=rt,
    at_risk_before=vapply(rt,function(t)sum(d$time_hours[ix]>=t),integer(1)))
 }
 point<-rbindlist(point);seed<-if(an=='PRIMARY_CARE_END_CENSORED')2026072906L else 2026072941L
 set.seed(seed);bt<-matrix(NA_real_,1000L,nrow(point))
 if(an=='PRIMARY_CARE_END_CENSORED')gb<-array(NA_real_,dim=c(1000L,length(grid),length(groups)))
 for(b in 1:1000){freq<-tabulate(sample.int(length(uid),length(uid),replace=TRUE),nbins=length(uid));weights<-freq[pid];k<-0L
  for(gi in seq_along(groups)){g<-groups[gi];ix<-if(g=='Overall')seq_len(nrow(d)) else which(d$dose_group==g)
   stopifnot(sum(weights[ix])>0);z<-aj(d$time_hours[ix],st[ix],weights[ix])
   for(cause in causes){bt[b,k+1:3]<-stepat(z,horizons,cause);k<-k+3L}
   if(an=='PRIMARY_CARE_END_CENSORED')gb[b,,gi]<-stepat(z,grid,'restart')
  }
  if(b%%250L==0){cat(an,'BOOTSTRAP',b,'/1000; elapsed seconds',round(as.numeric(difftime(Sys.time(),started,units='secs')),1),'\n');flush.console()}
 }
 stopifnot(all(is.finite(bt)))
 point[,`:=`(lower95=apply(bt,2,quantile,.025),upper95=apply(bt,2,quantile,.975),bootstrap_replicates=1000L,seed=seed)]
 estimates[[length(estimates)+1L]]<-point
 boots[[length(boots)+1L]]<-data.table(analysis=an,bootstrap=rep(1:1000,nrow(point)),table_row=rep(seq_len(nrow(point)),each=1000),value=as.vector(bt))
 fwrite(point[,.(table_row=.I,dose_group,cause,time_hours)],file.path(a,paste0(if(an=='PRIMARY_CARE_END_CENSORED')'primary' else 's20','_bootstrap_column_key.csv')))
 if(an=='PRIMARY_CARE_END_CENSORED'){
  stopifnot(all(is.finite(gb)))
  for(gi in seq_along(groups))for(j in seq_along(horizons))band_delta<-max(band_delta,max(abs(bt[,(gi-1L)*9L+j]-gb[,match(horizons[j],grid),gi])))
  stopifnot(band_delta<1e-14)
  fullband<-rbindlist(lapply(seq_along(groups),function(gi){g<-groups[gi];ix<-if(g=='Overall')seq_len(nrow(d)) else which(d$dose_group==g)
   data.table(cohort='MIMIC',version='minimal_1958',analysis=an,dose_group=g,time_hours=grid,
    estimate=stepat(aj(d$time_hours[ix],st[ix]),grid,'restart'),lower95=apply(gb[,,gi],2,quantile,.025),
    median_bootstrap=apply(gb[,,gi],2,median),upper95=apply(gb[,,gi],2,quantile,.975),bootstrap_replicates=1000L,seed=seed)}))
  fullboot<-rbindlist(lapply(seq_along(groups),function(gi)data.table(dose_group=groups[gi],bootstrap=rep(1:1000,length(grid)),time_hours=rep(grid,each=1000),restart_CIF=as.vector(gb[,,gi]))))
 }
}
estimates<-rbindlist(estimates);curves<-rbindlist(curves);checks<-rbindlist(checks);risks<-rbindlist(risks)
comparison<-merge(fullband[time_hours%in%horizons],estimates[analysis=='PRIMARY_CARE_END_CENSORED'&cause=='restart'],
 by=c('cohort','version','analysis','dose_group','time_hours'),suffixes=c('_grid','_horizon'))
stopifnot(nrow(comparison)==12L,max(abs(comparison$estimate_grid-comparison$estimate_horizon))<1e-14,
          max(abs(comparison$lower95_grid-comparison$lower95_horizon))<1e-14,max(abs(comparison$upper95_grid-comparison$upper95_horizon))<1e-14,
          all(estimates$lower95>=0),all(estimates$upper95<=1),all(curves$time_hours>=4),all(risks$time_hours>=4))
fwrite(estimates,file.path(a,'cumulative_incidence_12_24_72h.csv'))
fwrite(curves,file.path(a,'cumulative_incidence_curves.csv'))
fwrite(risks,file.path(a,'risk_sets.csv'))
fwrite(rbindlist(boots),file.path(a,'patient_cluster_bootstrap.csv'))
fwrite(checks,file.path(a,'AJ_independent_checks.csv'))
fwrite(fullband,file.path(a,'mimic_primary_full_grid_CI.csv'))
fwrite(fullboot,file.path(a,'mimic_primary_full_grid_bootstrap.csv'))
verify_sources();stopifnot(digest(file=input,algo='sha256')==input_hash)
qc<-list(PASS=TRUE,events=nrow(d),patients=uniqueN(d$subject_id),bootstrap_success_each=1000L,bootstrap_analyses=2L,
 primary_seed=2026072906L,S20_seed=2026072941L,workers=1L,data_table_threads=getDTthreads(),
 reference_checks=nrow(checks),reference_max_difference=max(checks$max_absolute_difference),
 entry_hours=4L,horizon_hours=72L,before_q_curve_rows=0L,basic_tests=length(tests),
 immediate_restart_at_q=sum(d$status_main==1&d$time_hours==4),
 full_grid_start=4,full_grid_end=72,full_grid_step=.25,full_grid_points=length(grid),
 full_grid_same_draw_horizon_max_difference=band_delta,
 confidence_band='pointwise patient-cluster percentile; not simultaneous; not three-point interpolation',
 source_hashes_unchanged=TRUE,input_hash_unchanged=TRUE,input_sha256=input_hash,
 external_cohort_analyses=0L,historical_cohort_analyses=0L,new_regression_models=0L,event_readjudications=0L,
 elapsed_seconds=as.numeric(difftime(Sys.time(),started,units='secs')))
write_json(qc,file.path(a,'cif_qc.json'),pretty=TRUE,auto_unbox=TRUE,digits=16)
capture.output(sessionInfo(),file=file.path(a,'R_session_info.txt'))
print(estimates[cause=='restart']);print(qc)
