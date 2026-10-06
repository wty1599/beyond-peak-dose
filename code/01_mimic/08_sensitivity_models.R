options(stringsAsFactors=FALSE,warn=1,digits=17)
Sys.setenv(OMP_NUM_THREADS=1,OPENBLAS_NUM_THREADS=1,MKL_NUM_THREADS=1)
suppressPackageStartupMessages({library(data.table);library(Hmisc);library(jsonlite);library(digest)})
setDTthreads(1L)
source('code/00_setup/paths.R')
out<-private_output('mimic_sensitivity')
dir.create(file.path(out,'restricted_work'),showWarnings=FALSE)
started<-Sys.time();seed<-202609222L
logline<-function(...) {cat(format(Sys.time(),'%Y-%m-%d %H:%M:%S'),paste(...,collapse=' '),'\n');flush.console()}
status<-function(state,...)write_json(c(list(state=state,time=as.character(Sys.time())),list(...)),file.path(out,'RUN_STATUS.json'),pretty=TRUE,auto_unbox=TRUE,digits=17)
status('INPUT_PRECHECK')
paths<-c(cohort=release_path('mimic_selected_cohort'),model_input=release_path('mimic_model_input'),raw_input=release_path('mimic_observed_input'),knots=release_path('frozen_knots'),winsor=release_path('frozen_winsor'),clinical_input=release_path('mimic_clinical_input'),t2_script='code/01_mimic/03_association_models.R',t2_contrasts=file.path(private_output('mimic_association'),'standardized_contrasts_with_ci.csv'))
hashes<-vapply(paths,function(p)digest(file=p,algo='sha256'),character(1))
vars<-c('predictor_mechanical_ventilation_t0','predictor_rrt_t0','predictor_age','predictor_female','predictor_episode_peak_nee','predictor_positive_nee_hours','predictor_sofa_24h')
spl<-vars[c(3,5,6,7)]
all_dat<-fread(paths['model_input'],select=c('stay_id','subject_id','outcome_restart_72h',vars))
raw<-fread(paths['raw_input'],select=c('stay_id',vars))
cohort<-fread(paths['cohort'],select=c('stay_id','subject_id','limitation_state','restart_72h','death_72h','composite_72h'))
clin<-fread(paths['clinical_input'])
stopifnot(nrow(all_dat)==1958L,uniqueN(all_dat$subject_id)==1856L,uniqueN(all_dat$stay_id)==1958L,
 nrow(cohort)==1958L,uniqueN(cohort$stay_id)==1958L,nrow(raw)==1958L,uniqueN(raw$stay_id)==1958L)
ci<-match(all_dat$stay_id,cohort$stay_id);ri<-match(all_dat$stay_id,raw$stay_id)
stopifnot(!anyNA(ci),!anyNA(ri),all(all_dat$subject_id==cohort$subject_id[ci]),all(all_dat$outcome_restart_72h==cohort$restart_72h[ci]),
 all(cohort$limitation_state%in%c('LIMITATION_NONE','LIMITATION_POST','LIMITATION_UNKNOWN')),sum(all_dat$outcome_restart_72h)==680L)
all_dat[,`:=`(limitation_state=sub('^LIMITATION_','',cohort$limitation_state[ci]),composite_72h=cohort$composite_72h[ci])]
fwrite(data.table(source_value=c('LIMITATION_NONE','LIMITATION_POST','LIMITATION_UNKNOWN'),analysis_label=c('NONE','POST','UNKNOWN'),operation='remove fixed display prefix only; no reclassification'),file.path(out,'limitation_label_mapping.csv'))
expected_composite<-ifelse(cohort$restart_72h==1 | (!is.na(cohort$death_72h)&cohort$death_72h==1),1L,ifelse(cohort$restart_72h==0 & !is.na(cohort$death_72h)&cohort$death_72h==0,0L,NA_integer_))
stopifnot(identical(is.na(expected_composite),is.na(cohort$composite_72h)),all(expected_composite[!is.na(expected_composite)]==cohort$composite_72h[!is.na(expected_composite)]),
 sum(!is.na(cohort$composite_72h))==1945L,sum(cohort$composite_72h,na.rm=TRUE)==713L)
qc<-data.table(variable=vars,missing=vapply(all_dat[,..vars],function(x)sum(is.na(x)),integer(1)),nonfinite=vapply(all_dat[,..vars],function(x)sum(!is.finite(x)),integer(1)))
fwrite(qc,file.path(out,'input_missingness.csv'));stopifnot(all(qc$missing==0),all(qc$nonfinite==0))
kt<-fread(paths['knots']);wt<-fread(paths['winsor'])
knots<-lapply(spl,function(v)kt[predictor==v][order(knot_index)]$knot_value);names(knots)<-spl
for(v in spl)stopifnot(identical(pmin(pmax(raw[[v]][ri],wt[predictor==v]$lower_p005),wt[predictor==v]$upper_p995),all_dat[[v]]))
fwrite(kt[predictor%in%spl],file.path(out,'frozen_knots_used.csv'));fwrite(wt[predictor%in%spl],file.path(out,'frozen_winsor_used.csv'))
lim<-all_dat[,.(events=.N,patients=uniqueN(subject_id),restarts=sum(outcome_restart_72h),composite_evaluable=sum(!is.na(composite_72h)),composite_positive=sum(composite_72h,na.rm=TRUE),composite_unknown=sum(is.na(composite_72h))),by=limitation_state]
fwrite(lim,file.path(out,'treatment_limitation_counts.csv'))
stopifnot(nrow(clin)==1958L,uniqueN(clin$stay_id)==1958L,sum(clin$state=='RESTART')==680L)
cj<-match(all_dat$stay_id,clin$stay_id);stopifnot(!anyNA(cj),all((clin$state[cj]=='RESTART')==(all_dat$outcome_restart_72h==1)))
restarted<-clin[state=='RESTART'];tt<-restarted$time_hours
stopifnot(length(tt)==680L,all(is.finite(tt)),all(tt>=4 & tt<=72))
restart_summary<-data.table(population='MIMIC minimal correction; observed first restart events only',events=length(tt),patients=uniqueN(restarted$subject_id),time_origin='original t0',unit='hours',q1=as.numeric(quantile(tt,.25,type=7)),median=median(tt),q3=as.numeric(quantile(tt,.75,type=7)),quantile_type=7L)
fwrite(restart_summary,file.path(out,'restart_time_among_restarted.csv'))
fnames<-c('make_x','lp_check','fit_model','summarize_contrasts')
expressions<-parse(paths['t2_script']);got<-character()
for(e in expressions)if(is.call(e)&&identical(e[[1]],as.name('<-'))&&is.symbol(e[[2]])&&as.character(e[[2]])%in%fnames){eval(e,envir=.GlobalEnv);got<-c(got,as.character(e[[2]]))}
stopifnot(setequal(got,fnames))
fwrite(data.table(function_name=fnames,sha256=vapply(fnames,function(n)digest(paste(deparse(get(n)),collapse='\n'),algo='sha256',serialize=FALSE),character(1)),source_path=paths['t2_script'],operation='function definition reused unchanged'),file.path(out,'T2_function_reuse.csv'))
contrasts<-data.table(contrast=c('MV_1_vs_0','RRT_1_vs_0','duration_48_vs_24_h','peak_0.30_vs_0.15'),variable=vars[c(1,2,6,5)],low=c(0,0,24,.15),high=c(1,1,48,.30))
selection<-list(restart_exclude_POST=which(all_dat$limitation_state!='POST'),restart_known_status=which(all_dat$limitation_state!='UNKNOWN'),composite_evaluable=which(!is.na(all_dat$composite_72h)))
endpoints<-c(restart_exclude_POST='outcome_restart_72h',restart_known_status='outcome_restart_72h',composite_evaluable='composite_72h')
expected_n<-c(restart_exclude_POST=1437L,restart_known_status=1324L,composite_evaluable=1945L)
expected_patients<-c(restart_exclude_POST=1368L,restart_known_status=1263L,composite_evaluable=1843L)
expected_positive<-c(restart_exclude_POST=415L,restart_known_status=489L,composite_evaluable=713L)
maincoef<-maincon<-mainqc<-influence<-seps<-supports<-samples<-bootcoefs<-bootcons<-bootqcs<-models<-list()
failures<-data.table(replicate=integer(),model=character(),reason=character())
fwrite(failures,file.path(out,'bootstrap_failures.csv'))
for(m in names(selection)){
 dat<-all_dat[selection[[m]]];y<-dat[[endpoints[m]]];x<-make_x(dat,FALSE)
 stopifnot(ncol(x)==13L,all(y%in%0:1),!anyNA(y),nrow(dat)==expected_n[m],uniqueN(dat$subject_id)==expected_patients[m],sum(y)==expected_positive[m])
 samples[[m]]<-data.table(model=m,outcome=endpoints[m],events=nrow(dat),patients=uniqueN(dat$subject_id),positive=sum(y),negative=sum(y==0),excluded_from1958=1958L-nrow(dat),NONE=sum(dat$limitation_state=='NONE'),POST=sum(dat$limitation_state=='POST'),UNKNOWN=sum(dat$limitation_state=='UNKNOWN'),standardization_target='fixed own analysis subset',seed=seed)
 fwrite(rbindlist(samples),file.path(out,'analysis_populations.csv'))
 sp<-rbindlist(c(lapply(vars[1:6],function(v){ind<-c(1L,grep(paste0('^',v,'($|_rcs)'),colnames(x)));lp_check(x[,ind,drop=FALSE],v)}),list(lp_check(x,'full_core12'))))
 sp[,model:=m];seps[[m]]<-sp;fwrite(rbindlist(seps),file.path(out,'separation_linear_program.csv'))
 if(any(sp$solved!=1)||any(sp$separation)){status('STOP_SEPARATION',model=m);stop('Separation or inconclusive LP; author decision required')}
 su<-rbindlist(lapply(seq_len(nrow(contrasts)),function(j){v<-contrasts$variable[j];p<-quantile(dat[[v]],c(.05,.95),type=7);data.table(model=m,contrast=contrasts$contrast[j],value=c(contrasts$low[j],contrasts$high[j]),p5=p[1],p95=p[2],within_p5_p95=c(contrasts$low[j],contrasts$high[j])>=p[1]&c(contrasts$low[j],contrasts$high[j])<=p[2])}))
 supports[[m]]<-su;fwrite(rbindlist(supports),file.path(out,'contrast_support.csv'))
 if(any(!su$within_p5_p95)){status('STOP_SUPPORT',model=m);stop('Unexpected support failure; author decision required')}
 contrast_x<-setNames(list(lapply(seq_len(nrow(contrasts)),function(j){a<-copy(dat);b<-copy(dat);set(a,j=contrasts$variable[j],value=rep(contrasts$high[j],nrow(a)));set(b,j=contrasts$variable[j],value=rep(contrasts$low[j],nrow(b)));list(high=make_x(a,FALSE),low=make_x(b,FALSE))})),m)
 f<-fit_model(x,y)
 mainqc[[m]]<-data.table(model=m,events=length(y),patients=uniqueN(dat$subject_id),positive=sum(y),df=12L,rank=f$fit$rank,converged=f$fit$converged,iterations=f$fit$iter,min_prediction=min(f$fit$fitted.values),max_prediction=max(f$fit$fitted.values),near_separation=f$near,warnings=f$warnings)
 fwrite(rbindlist(mainqc),file.path(out,'main_fit_checks.csv'))
 if(f$bad||f$near){status('STOP_MAIN_FIT',model=m);stop('Main fit failure or near separation')}
 models[[m]]<-f$fit;maincoef[[m]]<-data.table(model=m,term=names(f$fit$coefficients),coefficient=as.numeric(f$fit$coefficients));maincon[[m]]<-summarize_contrasts(f$fit$coefficients,m)
 influence[[m]]<-data.table(model=m,rank=1:10,cooks_distance=sort(cooks.distance(f$fit),decreasing=TRUE)[1:10])
 fwrite(rbindlist(maincoef),file.path(out,'coefficients_point.csv'));fwrite(rbindlist(maincon),file.path(out,'standardized_contrasts_point.csv'))
 groups<-split(seq_len(nrow(dat)),dat$subject_id);RNGkind('Mersenne-Twister','Inversion','Rejection');set.seed(seed)
 bc<-bt<-bq<-list();nf<-0L
 for(b in 1:1000){
  ix<-unlist(groups[sample.int(length(groups),length(groups),replace=TRUE)],use.names=FALSE)
  ff<-tryCatch(fit_model(x[ix,,drop=FALSE],y[ix]),error=function(e)list(error=conditionMessage(e)))
  if(!is.null(ff$error)){bad<-TRUE;near<-FALSE;why<-ff$error;conv<-FALSE;rk<-NA_integer_;it<-NA_integer_;lo<-hi<-NA_real_}else{bad<-ff$bad;near<-ff$near;why<-ff$warnings;conv<-ff$fit$converged;rk<-ff$fit$rank;it<-ff$fit$iter;lo<-min(ff$fit$fitted.values);hi<-max(ff$fit$fitted.values)}
  bq[[b]]<-data.table(model=m,replicate=b,success=!bad&&!near,converged=conv,rank=rk,iterations=it,bootstrap_events=length(ix),positive=sum(y[ix]),min_prediction=lo,max_prediction=hi,near_separation=near,warnings=why)
  if(bad||near){
   nf<-nf+1L;failures<-rbind(failures,data.table(replicate=b,model=m,reason=if(near)paste('Near separation',why)else why));fwrite(failures,file.path(out,'bootstrap_failures.csv'))
   if(near||nf>10L){fwrite(rbindlist(c(bootqcs,list(rbindlist(bq)))),file.path(out,'bootstrap_fit_checks.csv'));status('STOP_BOOTSTRAP',model=m,replicate=b,failures=nf);stop('Bootstrap failure gate')}
  }else{
   bc[[b]]<-data.table(model=m,replicate=b,term=names(ff$fit$coefficients),coefficient=as.numeric(ff$fit$coefficients))
   z<-summarize_contrasts(ff$fit$coefficients,m);z[,replicate:=b];bt[[b]]<-z
  }
  if(b%%100==0){logline(m,b,'/1000; failures',nf);status('BOOTSTRAP_RUNNING',model=m,completed=b,failures=nf)}
 }
 bootcoefs[[m]]<-rbindlist(bc);bootcons[[m]]<-rbindlist(bt);bootqcs[[m]]<-rbindlist(bq)
 fwrite(rbindlist(bootcoefs),file.path(out,'bootstrap_coefficients.csv'));fwrite(rbindlist(bootcons),file.path(out,'bootstrap_standardized_contrasts.csv'));fwrite(rbindlist(bootqcs),file.path(out,'bootstrap_fit_checks.csv'))
 logline(m,'completed')
}
fwrite(rbindlist(influence),file.path(out,'cooks_top10_no_identifiers.csv'))
saveRDS(list(models=models,knots=knots,contrasts=contrasts,seed=seed,targets='fixed respective analytic subset'),file.path(out,'restricted_work/core_models.rds'))
bc<-rbindlist(bootcoefs);bt<-rbindlist(bootcons);bq<-rbindlist(bootqcs)
q<-function(x)quantile(x,c(.025,.975),type=7,names=FALSE)
cc<-bc[,.(coefficient_lower=q(coefficient)[1],coefficient_upper=q(coefficient)[2],successful_bootstrap=.N),by=.(model,term)]
cc<-merge(rbindlist(maincoef),cc,by=c('model','term'),sort=FALSE);fwrite(cc,file.path(out,'coefficients_with_bootstrap_ci.csv'))
long<-melt(bt,id.vars=c('model','contrast','replicate'),measure.vars=c('standardized_risk_low','standardized_risk_high','risk_difference_pp','adjusted_OR'),variable.name='measure')
cis<-long[,.(lower=q(value)[1],upper=q(value)[2],successful_bootstrap=.N),by=.(model,contrast,measure)]
pts<-melt(rbindlist(maincon),id.vars=c('model','contrast'),variable.name='measure',value.name='estimate')
fin<-merge(pts,cis,by=c('model','contrast','measure'),sort=FALSE);fwrite(fin,file.path(out,'standardized_contrasts_with_ci.csv'))
base<-fread(paths['t2_contrasts'])[model=='core12'];base[,`:=`(model='original_core12_restart',source='existing T2, directly reused')]
comparison<-rbind(base,copy(fin)[,source:='current requested supplementary analysis'],fill=TRUE)
fwrite(comparison,file.path(out,'core_models_comparison_long.csv'))
after<-vapply(paths,function(p)digest(file=p,algo='sha256'),character(1))
fwrite(data.table(role=names(paths),path=unname(paths),sha256_before=unname(hashes),sha256_after=unname(after),unchanged=unname(hashes==after)),file.path(out,'source_hash_comparison.csv'))
stopifnot(all(after==hashes),nrow(fin)==48L,nrow(bq)==3000L)
capture.output(sessionInfo(),file=file.path(out,'R_session_info.txt'))
summary<-list(computation_complete=TRUE,cohort_events=1958L,cohort_patients=1856L,restarts=680L,
 models=3L,df_each=12L,seed=seed,bootstrap_attempts_each=1000L,bootstrap_fits=nrow(bq),failures=nrow(failures),
 source_hashes_unchanged=TRUE,protected_after='PENDING',independent_audit='PENDING',model_refits_27df=0L,
 event_readjudications=0L,external_analyses=0L,old_manuscripts_modified=FALSE,elapsed_seconds=as.numeric(difftime(Sys.time(),started,units='secs')))
write_json(summary,file.path(out,'COMPUTATION_QC.json'),pretty=TRUE,auto_unbox=TRUE,digits=17)
status('COMPLETED_PENDING_AUDIT',bootstrap_model_fits=nrow(bq),failures=nrow(failures))
print(rbindlist(samples));print(restart_summary);print(fin[measure%in%c('risk_difference_pp','adjusted_OR')]);print(summary)
