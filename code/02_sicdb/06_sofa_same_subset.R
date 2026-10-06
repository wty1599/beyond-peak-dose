options(stringsAsFactors=FALSE, warn=1)
Sys.setenv(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
suppressPackageStartupMessages({library(data.table);library(Hmisc);library(jsonlite);library(digest)})
setDTthreads(1L)
source('code/00_setup/paths.R')
root<-private_output('sicdb_sofa_subset');a<-root;w<-file.path(root,'private');dir.create(w,showWarnings=FALSE)
started<-Sys.time()
kpath<-release_path('frozen_knots');wpath<-release_path('frozen_winsor')
kt <- fread(kpath);knots <- split(kt$knot_value,kt$predictor);win <- fread(wpath)
base <- c('predictor_age','predictor_female','predictor_episode_peak_nee','predictor_positive_nee_hours','predictor_mechanical_ventilation_t0','predictor_rrt_t0')
sf <- 'predictor_sofa_24h'
readinput <- function(fn,n,patients,events,sofa,prepared){
  req <- c(base,if(sofa)sf)
  selected <- c('CaseID','PatientID','outcome_restart_72h',req,'recorded_crrt_active_pre6h','predictor_rrt_t0_recorded')
  z <- fread(release_path(if(sofa)'sicdb_sofa_input' else 'sicdb_association_input'),select=selected)
  stopifnot(nrow(z)==n,uniqueN(z$PatientID)==patients,sum(z$outcome_restart_72h)==events,!anyDuplicated(z$CaseID))
  stopifnot(all(z$recorded_crrt_active_pre6h==z$predictor_rrt_t0_recorded),!anyNA(z$recorded_crrt_active_pre6h))
  z[,predictor_rrt_t0:=recorded_crrt_active_pre6h]
  stopifnot(all(vapply(z[,..req],function(x)all(is.finite(x)),logical(1))),!anyNA(z[,..req]))
  for(k in c('predictor_female','predictor_mechanical_ventilation_t0','predictor_rrt_t0','outcome_restart_72h'))stopifnot(all(z[[k]]%in%c(0,1)))

  z
}
raw1553 <- readinput('S30_noSOFA12df_inputs.csv',1553L,1523L,517L,FALSE,'M12_S30/restricted_work/S30_association_inputs.csv')
raw846 <- readinput('SOFA15df_inputs.csv',846L,829L,315L,TRUE,'M15_SOFA/restricted_work/organ_support_association_inputs.csv')
stopifnot(all(raw846$CaseID%in%raw1553$CaseID))
fwrite(rbindlist(lapply(c('S30_1553','SOFA_subset_846'),function(nm){z<-if(nm=='S30_1553')raw1553 else raw846;data.table(sample=nm,variable=names(z)[startsWith(names(z),'predictor_')],missing=vapply(z[,startsWith(names(z),'predictor_'),with=FALSE],function(x)sum(is.na(x)),integer(1)))})),file.path(a,'required_predictor_missingness.csv'))
samples <- rbindlist(lapply(c('S30_1553','SOFA_subset_846'),function(nm){z<-if(nm=='S30_1553')raw1553 else raw846;qq<-quantile(z$predictor_positive_nee_hours,c(.25,.5,.75));data.table(sample=nm,n=nrow(z),patients=uniqueN(z$PatientID),restarts=sum(z$outcome_restart_72h),duration_q1=qq[1],duration_median=qq[2],duration_q3=qq[3],duration_IQR=qq[3]-qq[1],events_at_least_400=sum(z$outcome_restart_72h)>=400,duration_IQR_at_least_30h=unname(qq[3]-qq[1])>=30,MV_positive=sum(z$predictor_mechanical_ventilation_t0),recorded_CRRT_positive=sum(z$predictor_rrt_t0))}))
fwrite(samples,file.path(a,'sample_descriptives_and_gates.csv'))
exposures<-c('predictor_mechanical_ventilation_t0','predictor_rrt_t0','predictor_positive_nee_hours','predictor_episode_peak_nee')
values<-list(c(0,1),c(0,1),c(24,48),c(.15,.30))
supports<-rbindlist(lapply(c('S30_1553','SOFA_subset_846'),function(nm){z<-if(nm=='S30_1553')raw1553 else raw846;rbindlist(lapply(1:4,function(i){xx<-z[[exposures[i]]];qq<-quantile(xx,c(.05,.95));data.table(sample=nm,variable=exposures[i],value=values[[i]],p05=qq[1],p95=qq[2],percent_strictly_below=vapply(values[[i]],function(t)100*mean(xx<t),numeric(1)),percent_at_or_below=vapply(values[[i]],function(t)100*mean(xx<=t),numeric(1)),within_p05_p95=values[[i]]>=qq[1]&values[[i]]<=qq[2],scale=if(i<3)'binary_state_not_continuous_support_test' else 'continuous')}))}))
fwrite(supports,file.path(a,'contrast_support.csv'))
if(any(!supports[scale=='continuous']$within_p05_p95))stop('STOP: continuous contrast outside P5-P95; author adjudication required')
variablecheck<-rbindlist(lapply(c('S30_1553','SOFA_subset_846'),function(nm){z<-if(nm=='S30_1553')raw1553 else raw846;req<-c(base,if(nm=='SOFA_subset_846')sf);rbindlist(lapply(req,function(k)rbindlist(lapply(0:1,function(y){xx<-z[outcome_restart_72h==y][[k]];data.table(sample=nm,variable=k,outcome=y,n=length(xx),minimum=min(xx),maximum=max(xx),binary_zero=if(all(xx%in%c(0,1)))sum(xx==0) else NA_integer_,binary_one=if(all(xx%in%c(0,1)))sum(xx==1) else NA_integer_)}))))}))
fwrite(variablecheck,file.path(a,'per_variable_separation_screen.csv'))
clip <- function(k,x){ii<-match(k,win$predictor);if(is.na(ii))x else pmin(pmax(x,win$lower_p005[ii]),win$upper_p995[ii])}
prepare <- function(z){q<-copy(z);for(k in c(base,sf))if(k%in%names(q))set(q,j=k,value=clip(k,q[[k]]));q}
v1553<-prepare(raw1553);v<-prepare(raw846)
basis<-function(k,z){m<-Hmisc::rcspline.eval(z,knots=knots[[k]],inclx=TRUE,norm=2);colnames(m)<-paste0(k,'_rcs',1:3);m}
design<-function(z,sofa=FALSE){xx<-cbind('(Intercept)'=1,basis('predictor_age',z$predictor_age),predictor_female=z$predictor_female,basis('predictor_episode_peak_nee',z$predictor_episode_peak_nee),basis('predictor_positive_nee_hours',z$predictor_positive_nee_hours));if(sofa)xx<-cbind(xx,basis(sf,z[[sf]]));cbind(xx,predictor_mechanical_ventilation_t0=z$predictor_mechanical_ventilation_t0,predictor_rrt_t0=z$predictor_rrt_t0)}
xlist<-list(core12_same846=design(v),sofa15_same846=design(v,TRUE))
stopifnot(ncol(xlist[[1]])==13L,ncol(xlist[[2]])==16L)
source('code/02_sicdb/separation_check.R')
lpchecks<-rbindlist(lapply(c('S30_1553','core12_same846','sofa15_same846'),function(nm){xx<-if(nm=='S30_1553')design(v1553) else xlist[[nm]];yy<-if(nm=='S30_1553')v1553$outcome_restart_72h else v$outcome_restart_72h;separation_check(xx,yy,nm)}))
fwrite(lpchecks,file.path(a,'separation_linear_program_check.csv'))
stopifnot(all(lpchecks$pass))
fitone<-function(xx,yy){warnings_seen<-character();f<-withCallingHandlers(glm.fit(xx,yy,family=binomial(),control=glm.control(maxit=100,epsilon=1e-8)),warning=function(e){warnings_seen<<-c(warnings_seen,conditionMessage(e));invokeRestart('muffleWarning')});
  if(!f$converged||f$rank!=ncol(xx)||any(!is.finite(f$coefficients)))stop('Nonconvergence/rank/nonfinite coefficients')
  if(any(grepl('0 or 1|separat',warnings_seen))||min(f$fitted.values)<1e-8||max(f$fitted.values)>1-1e-8)stop('Potential separation or near separation')
  f$warnings_seen<-warnings_seen;f}
standardize<-function(xx,beta){lp<-drop(xx%*%beta);rbindlist(lapply(1:4,function(i){jj<-which(startsWith(colnames(xx),exposures[i]));nb<-if(i<3)matrix(values[[i]],ncol=1) else basis(exposures[i],clip(exposures[i],values[[i]]));base_lp<-lp-drop(xx[,jj,drop=FALSE]%*%beta[jj]);add<-drop(nb%*%beta[jj]);rr<-colMeans(plogis(outer(base_lp,add,'+')));data.table(variable=exposures[i],reference=values[[i]][1],comparison=values[[i]][2],risk_low=rr[1],risk_high=rr[2],risk_difference_pp=100*(rr[2]-rr[1]),OR=exp(add[2]-add[1]))}))}
yy<-v$outcome_restart_72h
fits<-lapply(xlist,fitone,yy=yy)
pointdiag<-rbindlist(lapply(names(fits),function(nm){f<-fits[[nm]];data.table(model=nm,n=length(yy),predictor_df=f$rank-1L,converged=f$converged,iterations=f$iter,min_probability=min(f$fitted.values),max_probability=max(f$fitted.values),warnings=paste(f$warnings_seen,collapse=' | '))}))
fwrite(pointdiag,file.path(a,'point_model_diagnostics.csv'))
pointct<-rbindlist(lapply(names(fits),function(nm){z<-standardize(xlist[[nm]],fits[[nm]]$coefficients);z[,model:=nm];z}))
pointcoef<-rbindlist(lapply(names(fits),function(nm)data.table(model=nm,term=names(fits[[nm]]$coefficients),estimate=unname(fits[[nm]]$coefficients))))
influence<-rbindlist(lapply(names(fits),function(nm){f<-fits[[nm]];h<-rowSums(qr.Q(f$qr)^2);pearson<-(yy-f$fitted.values)/sqrt(f$fitted.values*(1-f$fitted.values));cooks<-pearson^2*h/(f$rank*(1-h)^2);data.table(model=nm,descending_rank=1:10,cook_distance=sort(cooks,decreasing=TRUE)[1:10])}))
fwrite(influence,file.path(a,'top10_cook_magnitudes_no_identifiers.csv'))
qq<-quantile(raw846[[sf]],c(1/3,2/3));bins<-cut(raw846[[sf]],breaks=c(-Inf,qq,Inf),include.lowest=TRUE,right=TRUE)
rrtd<-data.table(sofa_tertile=as.character(bins),CRRT=raw846$predictor_rrt_t0)[,.(events=.N,recorded_CRRT_positive=sum(CRRT),recorded_CRRT_negative=sum(CRRT==0)),by=sofa_tertile]
rrtd[,`:=`(tertile_cut1=unname(qq[1]),tertile_cut2=unname(qq[2]))];fwrite(rrtd,file.path(a,'recorded_CRRT_by_SOFA_tertile.csv'))
ids<-split(seq_len(nrow(v)),v$PatientID);SEED<-202609095L;set.seed(SEED)
bcoef<-list();bct<-list();diag<-list();failed<-list();warns<-list();success<-0L
for(attempt in seq_len(1000L)){
  ii<-unlist(ids[sample.int(length(ids),length(ids),replace=TRUE)],use.names=FALSE)
  current<-tryCatch(lapply(xlist,function(xx)fitone(xx[ii,,drop=FALSE],yy[ii])),error=function(e)e)
  if(inherits(current,'error')){failed[[length(failed)+1L]]<-data.table(attempt=attempt,reason=conditionMessage(current));fwrite(rbindlist(failed),file.path(a,'bootstrap_failures.csv'));stop('STOP: bootstrap fitting issue retained; no replacement or repair authorized')}
  success<-success+1L
  bcoef[[success]]<-rbindlist(lapply(names(current),function(nm)data.table(attempt=attempt,model=nm,term=names(current[[nm]]$coefficients),estimate=unname(current[[nm]]$coefficients))))
  bct[[success]]<-rbindlist(lapply(names(current),function(nm){z<-standardize(xlist[[nm]][ii,,drop=FALSE],current[[nm]]$coefficients);z[,`:=`(attempt=attempt,model=nm)];z}))
  diag[[success]]<-rbindlist(lapply(names(current),function(nm){f<-current[[nm]];data.table(attempt=attempt,model=nm,converged=f$converged,rank=f$rank,iterations=f$iter,n_resampled_events=length(ii),min_probability=min(f$fitted.values),max_probability=max(f$fitted.values),warnings=paste(f$warnings_seen,collapse=' | '))}))
  for(nm in names(current))if(length(current[[nm]]$warnings_seen))warns[[length(warns)+1L]]<-data.table(attempt=attempt,model=nm,warning=current[[nm]]$warnings_seen)
  if(attempt%%100L==0L){cat('T3 paired bootstrap',attempt,'/1000\n');flush.console()}
}
bc<-rbindlist(bcoef);bt<-rbindlist(bct)
fwrite(bc,file.path(a,'paired_bootstrap_coefficients.csv'));fwrite(bt,file.path(a,'paired_bootstrap_contrasts.csv'));fwrite(rbindlist(diag),file.path(a,'bootstrap_convergence.csv'))
fwrite(if(length(failed))rbindlist(failed) else data.table(attempt=integer(),reason=character()),file.path(a,'bootstrap_failures.csv'))
fwrite(if(length(warns))rbindlist(warns) else data.table(attempt=integer(),model=character(),warning=character()),file.path(a,'bootstrap_warnings.csv'))
cc<-bc[,.(lower95=quantile(estimate,.025),upper95=quantile(estimate,.975)),by=.(model,term)]
pc<-merge(pointcoef,cc,by=c('model','term'),sort=FALSE);fwrite(pc,file.path(a,'same846_complete_coefficients_with_percentile_CI.csv'))
metrics<-c('risk_low','risk_high','risk_difference_pp','OR')
cts<-rbindlist(lapply(metrics,function(m){qq<-bt[,.(lower95=if(m=='OR')exp(quantile(log(get(m)),.025)) else quantile(get(m),.025),upper95=if(m=='OR')exp(quantile(log(get(m)),.975)) else quantile(get(m),.975)),by=.(model,variable,reference,comparison)];pp<-pointct[,.(model,variable,reference,comparison,estimate=get(m))];zz<-merge(pp,qq,by=c('model','variable','reference','comparison'),sort=FALSE);zz[,`:=`(measure=m,unit=if(m=='OR')'odds_ratio' else if(m=='risk_difference_pp')'percentage_points' else 'probability')];zz}))
fwrite(cts,file.path(a,'same846_standardized_contrasts.csv'))
saveRDS(list(raw846=raw846,prepared846=v,design=xlist,fits=fits,seed=SEED,bootstrap_attempts=1000L,knots=knots,winsor=win),file.path(w,'paired_same846_analysis.rds'))
capture.output(sessionInfo(),file=file.path(a,'sessionInfo.txt'))
