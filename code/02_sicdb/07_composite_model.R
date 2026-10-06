# Prepared without fitting. An explicit execution argument is required after S5.
options(stringsAsFactors=FALSE,warn=1,digits=17)
invisible(Sys.setlocale('LC_CTYPE','English_United States.utf8'))
Sys.setenv(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
suppressPackageStartupMessages({library(data.table);library(Hmisc);library(jsonlite);library(digest)})
setDTthreads(1L)
source('code/00_setup/paths.R')
out<-private_output('sicdb_composite')
started<-Sys.time();dir.create(file.path(out,'restricted_work'),showWarnings=FALSE)
SEED<-2026091030L;CURRENT_PHASE<-'INPUT_PREFLIGHT';CURRENT_ATTEMPT<-0L
src<-release_path('sicdb_clinical_input');oldsrc<-release_path('sicdb_association_input')
kpath<-release_path('frozen_knots');wpath<-release_path('frozen_winsor')
t3path<-'code/02_sicdb/06_sofa_same_subset.R';lppath<-'code/02_sicdb/separation_check.R'
sources<-c(src,oldsrc,kpath,wpath,t3path,lppath)
hashes<-function()data.table(path=sources,sha256=vapply(sources,function(f)digest(file=f,algo='sha256'),character(1)))
before<-hashes();fwrite(before,file.path(out,'S6_source_sha256_before.csv'))
capture.output(sessionInfo(),file=file.path(out,'S6_R_session_info.txt'))
write_status<-function(status,message='')write_json(list(status=status,phase=CURRENT_PHASE,attempt=CURRENT_ATTEMPT,message=message,seed=SEED,time=format(Sys.time(),tz='UTC',usetz=TRUE),elapsed_seconds=as.numeric(difftime(Sys.time(),started,units='secs'))),file.path(out,'S6_STATUS.json'),auto_unbox=TRUE,pretty=TRUE,digits=NA)
write_status('RUNNING')
failed<-list();warns<-list();bc<-list();bt<-list();bd<-list()
flush_boot<-function(){
 if(length(bc))fwrite(rbindlist(bc),file.path(out,'S6_bootstrap_coefficients.csv'))
 if(length(bt))fwrite(rbindlist(bt),file.path(out,'S6_bootstrap_contrasts.csv'))
 if(length(bd))fwrite(rbindlist(bd),file.path(out,'S6_bootstrap_diagnostics.csv'))
 fwrite(if(length(failed))rbindlist(failed)else data.table(phase=character(),attempt=integer(),reason=character()),file.path(out,'S6_failures.csv'))
 fwrite(if(length(warns))rbindlist(warns)else data.table(phase=character(),attempt=integer(),warning=character()),file.path(out,'S6_warnings.csv'))
}
load_definition<-function(path,name){
 exprs<-parse(path,encoding='UTF-8')
 chosen<-Filter(function(e)is.call(e)&&identical(e[[1]],as.name('<-'))&&is.symbol(e[[2]])&&as.character(e[[2]])==name,as.list(exprs))
 stopifnot(length(chosen)==1L,is.call(chosen[[1]][[3]]),identical(chosen[[1]][[3]][[1]],as.name('function')))
 eval(chosen[[1]],envir=.GlobalEnv)
 invisible(chosen[[1]])
}
tryCatch({
 d<-fread(src)
 stopifnot(nrow(d)==2044L,uniqueN(d$CaseID)==2044L,uniqueN(d$PatientID)==1996L)
 endpoints<-c('restart_binary','death_72h_binary','composite_restart_or_death')
 gate<-rbindlist(lapply(endpoints,function(k)data.table(endpoint=k,cases=nrow(d),known=sum(!is.na(d[[k]])),positive=sum(d[[k]]==1,na.rm=TRUE),unknown=sum(is.na(d[[k]])))))
 fwrite(gate,file.path(out,'S6_endpoint_gate_check.csv'))
 stopifnot(identical(as.integer(d$restart_binary),as.integer(d$outcome_restart_72h)))
 stopifnot(all(d$recorded_crrt_active_pre6h %in% c(0,1)),!anyNA(d$recorded_crrt_active_pre6h),
           identical(as.numeric(d$recorded_crrt_active_pre6h),as.numeric(d$predictor_rrt_t0_recorded)))
 d[,predictor_rrt_t0:=recorded_crrt_active_pre6h]
 req<-c('predictor_age','predictor_female','predictor_episode_peak_nee','predictor_positive_nee_hours','predictor_mechanical_ventilation_t0','predictor_rrt_t0')
 complete_required<-complete.cases(d[,..req])
 reproduced_old<-d[complete_required & !is.na(outcome_restart_72h)]
 old<-fread(oldsrc,select=c('CaseID','PatientID','outcome_restart_72h',req))
 stopifnot(nrow(reproduced_old)==1553L,uniqueN(reproduced_old$PatientID)==1523L,sum(reproduced_old$outcome_restart_72h)==517L,
           nrow(old)==1553L,!anyDuplicated(old$CaseID),setequal(reproduced_old$CaseID,old$CaseID))
 old_order<-match(reproduced_old$CaseID,old$CaseID)
 parity<-data.table(field=c('PatientID','outcome_restart_72h',req),identical_value=vapply(c('PatientID','outcome_restart_72h',req),function(k)isTRUE(all.equal(reproduced_old[[k]],old[[k]][old_order],check.attributes=FALSE)),logical(1)))
 fwrite(parity,file.path(out,'S6_original_S30_rule_parity.csv'));stopifnot(all(parity$identical_value))
 raw<-d[complete_required & !is.na(composite_restart_or_death)]
 stopifnot(nrow(raw)>0L,!anyDuplicated(raw$CaseID),all(raw$composite_restart_or_death %in% 0:1),
           all(vapply(raw[,..req],function(x)all(is.finite(x)),logical(1))))
 flow<-data.table(stage=c('current_all_cases','six_required_predictors_complete','original_restart_known_rule_reproduced','composite_known_core_sample','core_cases_with_restart_unknown','original_S30_cases_excluded_composite_unknown'),
                 cases=c(nrow(d),sum(complete_required),nrow(reproduced_old),nrow(raw),sum(is.na(raw$restart_binary)),sum(reproduced_old$CaseID %in% d[is.na(composite_restart_or_death)]$CaseID)))
 fwrite(flow,file.path(out,'S6_sample_flow.csv'))
 composite_positive<-raw[composite_restart_or_death==1]
 restarted_positive<-composite_positive[restart_binary==1]
 dead_positive<-composite_positive[death_72h_binary==1]
 stopifnot(!anyNA(restarted_positive$valid_restart_time),!anyNA(dead_positive$death))
 deaths_without_prior_restart<-sum(is.na(dead_positive$valid_restart_time)|dead_positive$valid_restart_time>=dead_positive$death)
 composition<-data.table(component=c('Composite-positive cases','Valid restart within 72 h','Death within 72 h without preceding valid restart'),
                         cases=c(nrow(composite_positive),nrow(restarted_positive),deaths_without_prior_restart),
                         denominator=nrow(raw),counting_rule=c('composite_restart_or_death == 1','composite == 1 and restart_binary == 1','composite == 1 and death_72h_binary == 1 and (valid_restart_time missing or valid_restart_time >= death)'))
 fwrite(composition,file.path(out,'S6_composite_component_counts.csv'))
 stopifnot(nrow(restarted_positive)+deaths_without_prior_restart==nrow(composite_positive))
 fwrite(data.table(variable=req,missing_in_all=vapply(d[,..req],function(x)sum(is.na(x)),integer(1)),missing_in_composite_sample=vapply(raw[,..req],function(x)sum(is.na(x)),integer(1))),file.path(out,'S6_missingness.csv'))
 for(k in c('predictor_female','predictor_mechanical_ventilation_t0','predictor_rrt_t0'))stopifnot(all(raw[[k]] %in% 0:1))
 kt<-fread(kpath);knots<-split(kt$knot_value,kt$predictor);win<-fread(wpath)
 clip<-function(k,z){j<-match(k,win$predictor);if(is.na(j))z else pmin(pmax(z,win$lower_p005[j]),win$upper_p995[j])}
 v<-copy(raw);cliprecords<-list()
 for(k in req){j<-match(k,win$predictor);if(!is.na(j)){cliprecords[[k]]<-data.table(variable=k,lower=win$lower_p005[j],upper=win$upper_p995[j],below_lower=sum(v[[k]]<win$lower_p005[j]),above_upper=sum(v[[k]]>win$upper_p995[j]));set(v,j=k,value=clip(k,v[[k]]))}}
 fwrite(rbindlist(cliprecords),file.path(out,'S6_frozen_winsor_application.csv'))
 fwrite(kt[predictor %in% req],file.path(out,'S6_frozen_knots_used.csv'))
 basis<-function(k,z){m<-Hmisc::rcspline.eval(z,knots=knots[[k]],inclx=TRUE,norm=2);colnames(m)<-paste0(k,'_rcs',1:3);m}
 x<-cbind('(Intercept)'=1,basis('predictor_age',v$predictor_age),predictor_female=v$predictor_female,
          basis('predictor_episode_peak_nee',v$predictor_episode_peak_nee),basis('predictor_positive_nee_hours',v$predictor_positive_nee_hours),
          predictor_mechanical_ventilation_t0=v$predictor_mechanical_ventilation_t0,predictor_rrt_t0=v$predictor_rrt_t0)
 y<-v$composite_restart_or_death
 stopifnot(ncol(x)==13L,qr(x)$rank==13L,!anyNA(x),all(is.finite(x)))
 fwrite(data.table(position=seq_len(ncol(x)),column=colnames(x)),file.path(out,'S6_design_columns.csv'))
 exposures<-c('predictor_mechanical_ventilation_t0','predictor_rrt_t0','predictor_positive_nee_hours','predictor_episode_peak_nee')
 values<-list(c(0,1),c(0,1),c(24,48),c(.15,.30))
 support<-rbindlist(lapply(3:4,function(i){xx<-raw[[exposures[i]]];qq<-quantile(xx,c(.05,.95),type=7);data.table(variable=exposures[i],value=values[[i]],p05=unname(qq[1]),p95=unname(qq[2]),percent_strictly_below=vapply(values[[i]],function(a)100*mean(xx<a),numeric(1)),percent_at_or_below=vapply(values[[i]],function(a)100*mean(xx<=a),numeric(1)),within_p05_p95=values[[i]]>=qq[1]&values[[i]]<=qq[2])}))
 fwrite(support,file.path(out,'S6_continuous_contrast_support.csv'))
 if(any(!support$within_p05_p95))stop('Continuous contrast support failure: author adjudication required')
 q<-quantile(raw$predictor_positive_nee_hours,c(.25,.5,.75),type=7,names=FALSE)
 counts<-data.table(cases=nrow(raw),patients=uniqueN(raw$PatientID),composite_positive=sum(y),composite_negative=sum(y==0),restart_unknown_included=sum(is.na(raw$restart_binary)),duration_q1=q[1],duration_median=q[2],duration_q3=q[3],duration_IQR=q[3]-q[1],events_at_least_400=sum(y)>=400,duration_IQR_at_least_30h=q[3]-q[1]>=30,MV_positive=sum(raw$predictor_mechanical_ventilation_t0),recorded_CRRT_positive=sum(raw$predictor_rrt_t0))
 fwrite(counts,file.path(out,'S6_sample_and_interpretation_gates.csv'))
 CURRENT_PHASE<-'SEPARATION_PREFLIGHT';write_status('RUNNING')
 load_definition(lppath,'separation_check')
 tests<-list(overlap=list(x=c(0,1,0,1),y=c(0,0,1,1),expected=FALSE),complete=list(x=c(-2,-1,1,2),y=c(0,0,1,1),expected=TRUE),quasi=list(x=c(-1,0,0,1),y=c(0,0,1,1),expected=TRUE))
 syn<-rbindlist(lapply(names(tests),function(nm){q<-tests[[nm]];z<-separation_check(cbind(1,q$x),q$y,nm);z[,expected:=q$expected];z}))
 fwrite(syn,file.path(out,'S6_separation_synthetic_tests.csv'));stopifnot(all(syn$LP_solved==1L),all(syn$separation_or_quasi_separation==syn$expected))
 lp<-rbindlist(c(list(separation_check(x,y,'composite_core12')),lapply(req,function(k){ii<-c(1L,which(startsWith(colnames(x),k)));separation_check(x[,ii,drop=FALSE],y,k)})))
 fwrite(lp,file.path(out,'S6_separation_diagnostics.csv'));stopifnot(all(lp$pass))
 load_definition(t3path,'fitone')
 standardize<-function(xx,beta){lp<-drop(xx%*%beta);rbindlist(lapply(1:4,function(i){jj<-which(startsWith(colnames(xx),exposures[i]));nb<-if(i<3)matrix(values[[i]],ncol=1) else basis(exposures[i],clip(exposures[i],values[[i]]));base_lp<-lp-drop(xx[,jj,drop=FALSE]%*%beta[jj]);add<-drop(nb%*%beta[jj]);rr<-colMeans(plogis(outer(base_lp,add,'+')));data.table(variable=exposures[i],reference=values[[i]][1],comparison=values[[i]][2],risk_low=rr[1],risk_high=rr[2],risk_difference_pp=100*(rr[2]-rr[1]),log_OR=add[2]-add[1],OR=exp(add[2]-add[1]))}))}
 diagone<-function(f,attempt,n,eventn)data.table(attempt=attempt,cases=n,composite_positive=eventn,rank=f$rank,predictor_df=f$rank-1L,converged=f$converged,iterations=f$iter,min_probability=min(f$fitted.values),max_probability=max(f$fitted.values),warnings=paste(f$warnings_seen,collapse=' | '))
 CURRENT_PHASE<-'MAIN_FIT';write_status('RUNNING')
 pointfit<-fitone(x,y);point<-standardize(x,pointfit$coefficients)
 fwrite(diagone(pointfit,0L,length(y),sum(y)),file.path(out,'S6_main_fit_diagnostics.csv'))
 if(length(pointfit$warnings_seen))warns[[1L]]<-data.table(phase='MAIN_FIT',attempt=0L,warning=pointfit$warnings_seen)
 pointcoef<-data.table(model='SICdb_composite_core12',term=names(pointfit$coefficients),estimate=unname(pointfit$coefficients))
 fwrite(pointcoef,file.path(out,'S6_coefficients_point.csv'));fwrite(point,file.path(out,'S6_contrasts_point.csv'))
 h<-rowSums(qr.Q(pointfit$qr)^2);pearson<-(y-pointfit$fitted.values)/sqrt(pointfit$fitted.values*(1-pointfit$fitted.values));cooks<-pearson^2*h/(pointfit$rank*(1-h)^2)
 fwrite(data.table(descending_rank=1:10,cook_distance=sort(cooks,decreasing=TRUE)[1:10]),file.path(out,'S6_influence_top10_no_identifiers.csv'))
 fwrite(raw[,c('CaseID','PatientID','restart_binary','composite_restart_or_death',req),with=FALSE],file.path(out,'restricted_work/S6_composite_association_inputs.csv'))
 ids<-split(seq_len(nrow(v)),v$PatientID);RNGkind('Mersenne-Twister','Inversion','Rejection');set.seed(SEED)
 CURRENT_PHASE<-'PATIENT_CLUSTER_BOOTSTRAP';write_status('RUNNING')
 for(attempt in seq_len(1000L)){
  CURRENT_ATTEMPT<-attempt
  ii<-unlist(ids[sample.int(length(ids),length(ids),replace=TRUE)],use.names=FALSE)
  f<-fitone(x[ii,,drop=FALSE],y[ii])
  bc[[attempt]]<-data.table(attempt=attempt,term=names(f$coefficients),estimate=unname(f$coefficients))
  ctmp<-standardize(x[ii,,drop=FALSE],f$coefficients);ctmp[,attempt:=attempt];bt[[attempt]]<-ctmp
  bd[[attempt]]<-diagone(f,attempt,length(ii),sum(y[ii]))
  if(length(f$warnings_seen))warns[[length(warns)+1L]]<-data.table(phase=CURRENT_PHASE,attempt=attempt,warning=f$warnings_seen)
  if(attempt%%100L==0L){flush_boot();write_status('RUNNING');cat('S6 bootstrap ',attempt,' / 1000\n');flush.console()}
 }
 CURRENT_PHASE<-'AGGREGATION';flush_boot()
 coefboot<-rbindlist(bc);contrastboot<-rbindlist(bt)
 ci<-coefboot[,.(lower95=quantile(estimate,.025,type=7,names=FALSE),upper95=quantile(estimate,.975,type=7,names=FALSE),successful_bootstrap=.N),by=term]
 coefout<-merge(pointcoef,ci,by='term',sort=FALSE);fwrite(coefout,file.path(out,'S6_coefficients.csv'))
 metrics<-c('risk_low','risk_high','risk_difference_pp','OR')
 result<-rbindlist(lapply(metrics,function(m){key<-if(m=='OR')'log_OR' else m;qq<-contrastboot[,.(lower95=quantile(get(key),.025,type=7,names=FALSE),upper95=quantile(get(key),.975,type=7,names=FALSE),successful_bootstrap=.N),by=.(variable,reference,comparison)];pp<-point[,.(variable,reference,comparison,estimate=get(m))];if(m=='OR')qq[,`:=`(lower95=exp(lower95),upper95=exp(upper95))];z<-merge(pp,qq,by=c('variable','reference','comparison'),sort=FALSE);z[,`:=`(model='SICdb_composite_core12',measure=m,unit=if(m=='OR')'odds_ratio'else if(m=='risk_difference_pp')'percentage_points'else'probability')];z}))
 fwrite(result,file.path(out,'S6_composite_core_contrasts.csv'))
 saveRDS(list(raw=raw[,c('CaseID','PatientID','restart_binary','composite_restart_or_death',req),with=FALSE],prepared=v[,c('CaseID','PatientID','composite_restart_or_death',req),with=FALSE],design=x,fit=pointfit,knots=knots,winsor=win,seed=SEED,exposures=exposures,values=values,standardization='empirical analysis sample; resampled event population within each patient bootstrap'),file.path(out,'restricted_work/S6_composite_core_artifact.rds'))
 after<-hashes();fwrite(after,file.path(out,'S6_source_sha256_after.csv'));check<-merge(before,after,by='path',suffixes=c('_before','_after'));check[,unchanged:=sha256_before==sha256_after];fwrite(check,file.path(out,'S6_source_hash_check.csv'));stopifnot(all(check$unchanged))
 capture.output(sessionInfo(),file=file.path(out,'S6_R_session_info.txt'))
 CURRENT_PHASE<-'COMPLETE'
 write_json(list(status='CALCULATION_COMPLETE_PENDING_INDEPENDENT_REVIEW',cases=nrow(raw),patients=uniqueN(raw$PatientID),composite_events=sum(y),predictor_df=12,design_columns=13,seed=SEED,bootstrap_attempts=1000,bootstrap_success=length(bt),bootstrap_failures=0,required_covariates_complete=TRUE,reproduced_original_S30_cases=1553,reproduced_original_S30_patients=1523,reproduced_original_S30_restarts=517,restart_known_not_required_for_new_sample=TRUE,restart_unknown_included=sum(is.na(raw$restart_binary)),CRRT_interface='recorded_crrt_active_pre6h; checked against predictor_rrt_t0_recorded',new_imputations=0,SOFA_included=FALSE,source_hashes_unchanged=TRUE,bootstrap_standardization='resampled event population',OR_interval='percentiles of log_OR then exponentiation',events_gate=counts$events_at_least_400,duration_IQR_gate=counts$duration_IQR_at_least_30h,post_hoc_supportive_composite_analysis=TRUE,elapsed_seconds=as.numeric(difftime(Sys.time(),started,units='secs'))),file.path(out,'S6_QC.json'),auto_unbox=TRUE,pretty=TRUE,digits=NA)
 write_status('CALCULATION_COMPLETE_PENDING_INDEPENDENT_REVIEW')
 cat('S6 COMPLETE, no AUC or block deletion analyses performed\n')
},error=function(e){
 failed[[length(failed)+1L]]<<-data.table(phase=CURRENT_PHASE,attempt=CURRENT_ATTEMPT,reason=conditionMessage(e));flush_boot()
 after<-hashes();fwrite(after,file.path(out,'S6_source_sha256_after.csv'));z<-merge(before,after,by='path',suffixes=c('_before','_after'));z[,unchanged:=sha256_before==sha256_after];fwrite(z,file.path(out,'S6_source_hash_check.csv'))
 write_status('STOP_AUTHOR_REVIEW_REQUIRED',conditionMessage(e));stop(e)
})
