options(stringsAsFactors=FALSE, warn=1, digits=17)
Sys.setenv(OMP_NUM_THREADS=1, OPENBLAS_NUM_THREADS=1, MKL_NUM_THREADS=1)
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(Hmisc))
suppressPackageStartupMessages(library(jsonlite))
setDTthreads(1)
source('code/00_setup/paths.R')
out <- private_output('mimic_association')
dir.create(file.path(out,'restricted_work'),recursive=TRUE,showWarnings=FALSE)
start <- Sys.time()
logline <- function(...) cat(format(Sys.time(),'%Y-%m-%d %H:%M:%S'),paste(...,collapse=' '),'\n')
status <- function(state,...) write_json(c(list(state=state,time=format(Sys.time(),tz='UTC',usetz=TRUE)),list(...)),file.path(out,'run_status.json'),auto_unbox=TRUE,pretty=TRUE,digits=NA)
status('STARTED')
seed <- 202609222L
logline('Seed',seed,'1000 patient-cluster bootstrap repetitions; fixed original 1958-event standardization target')
input <- release_path('mimic_model_input')
rawfile <- release_path('mimic_observed_input')
kfile <- release_path('frozen_knots')
wfile <- release_path('frozen_winsor')
vars <- c('predictor_mechanical_ventilation_t0','predictor_rrt_t0','predictor_age','predictor_female','predictor_episode_peak_nee','predictor_positive_nee_hours','predictor_sofa_24h')
spl <- vars[c(3,5,6,7)]
dat <- fread(input,select=c('subject_id','outcome_restart_72h',vars))
raw <- fread(rawfile,select=vars)
stopifnot(nrow(dat)==1958L,uniqueN(dat$subject_id)==1856L,sum(dat$outcome_restart_72h)==680L,all(dat$outcome_restart_72h %in% 0:1))
qc <- data.table(variable=vars,missing=vapply(dat[,..vars],function(x)sum(is.na(x)),integer(1)),nonfinite=vapply(dat[,..vars],function(x)sum(!is.finite(x)),integer(1)))
fwrite(qc,file.path(out,'input_missingness.csv'));stopifnot(all(qc$missing==0L),all(qc$nonfinite==0L))
for(v in vars[c(1,2,4)]) stopifnot(all(dat[[v]] %in% 0:1))
kt <- fread(kfile); wt <- fread(wfile)
knots <- lapply(spl,function(v)kt[predictor==v][order(knot_index)]$knot_value);names(knots)<-spl
for(v in spl) stopifnot(identical(pmin(pmax(raw[[v]],wt[predictor==v]$lower_p005),wt[predictor==v]$upper_p995),dat[[v]]))
fwrite(kt[predictor %in% spl],file.path(out,'frozen_knots_used.csv'))
fwrite(wt[predictor %in% spl],file.path(out,'frozen_winsor_used.csv'))
make_x <- function(d,extended=FALSE) {
  parts<-list('(Intercept)'=matrix(1,nrow(d),1,dimnames=list(NULL,'(Intercept)')))
  for(v in vars[1:6]) {
    b <- if(v %in% spl) Hmisc::rcspline.eval(d[[v]],knots=knots[[v]],inclx=TRUE,norm=2) else matrix(d[[v]],ncol=1)
    colnames(b) <- if(v %in% spl) paste0(v,'_rcs',1:3) else v
    parts[[v]]<-b
  }
  if(extended) {v<-vars[7];b<-Hmisc::rcspline.eval(d[[v]],knots=knots[[v]],inclx=TRUE,norm=2);colnames(b)<-paste0(v,'_rcs',1:3);parts[[v]]<-b}
  do.call(cbind,parts)
}
xs <- list(core12=make_x(dat),sofa15=make_x(dat,TRUE))
stopifnot(ncol(xs$core12)==13L,ncol(xs$sofa15)==16L)
y<-dat$outcome_restart_72h
# Linear programming: maximize summed signed margins, each margin >= 0,
# with L1 norm of the free separating direction bounded by 1.
lp_check <- function(x,label) {
  x<-sweep(x,2,apply(abs(x),2,max),'/');z<-x*(2*y-1);p<-ncol(z)
  a<-c(colSums(z),-colSums(z));zz<-cbind(z,-z)
  ans<-boot::simplex(a=a,A1=rbind(-zz,rep(1,2*p)),b1=c(rep(0,nrow(z)),1),maxi=TRUE,n.iter=50000,eps=1e-10)
  dir<-ans$soln[seq_len(p)]-ans$soln[p+seq_len(p)]
  mar<-drop(z %*% dir)
  data.table(block=label,solved=ans$solved,objective=ans$value,min_margin=min(mar),max_margin=max(mar),separation=as.logical(ans$solved==1 && ans$value>1e-7 && min(mar)> -1e-7))
}
sep <- list()
for(v in vars) {
  selected<-c(1L,grep(paste0('^',v,'($|_rcs)'),colnames(xs$sofa15)))
  sep[[v]]<-lp_check(xs$sofa15[,selected,drop=FALSE],v)
  logline('Separation block',v,'done')
}
for(m in names(xs)) {sep[[m]]<-lp_check(xs[[m]],m);logline('Separation full',m,'done')}
sep<-rbindlist(sep);fwrite(sep,file.path(out,'separation_linear_program.csv'))
if(any(sep$solved!=1)||any(sep$separation)) {status('STOP_SEPARATION',message='LP failure or separation requires author review');stop('LP failure or separation; no fit')}
covtab<-rbindlist(lapply(vars,function(v){rbindlist(lapply(0:1,function(a){xx<-dat[y==a][[v]];data.table(variable=v,outcome=a,n=length(xx),minimum=min(xx),q1=quantile(xx,.25),median=median(xx),q3=quantile(xx,.75),maximum=max(xx),zero=sum(xx==0),one=sum(xx==1))}))}))
fwrite(covtab,file.path(out,'variable_outcome_coverage.csv'))
cutpoints<-quantile(dat$predictor_sofa_24h,c(0,1/3,2/3,1),type=7)
stopifnot(length(unique(cutpoints))==4)
tertile<-cut(dat$predictor_sofa_24h,breaks=cutpoints,include.lowest=TRUE,right=TRUE)
rrt<-as.data.table(table(SOFA_tertile=tertile,RRT=dat$predictor_rrt_t0));setnames(rrt,'N','events')
rrt[,within_tertile_percent:=100*events/sum(events),by=SOFA_tertile]
rrt[,among_RRT_positive_percent:=ifelse(RRT=='1',100*events/sum(events[RRT=='1']),NA_real_)]
fwrite(rrt,file.path(out,'rrt_by_sofa_tertile.csv'))
fwrite(data.table(probability=c(0,1/3,2/3,1),SOFA_boundary=as.numeric(cutpoints)),file.path(out,'sofa_tertile_boundaries.csv'))
if(any(rrt$events==0)) {status('STOP_JOINT_SUPPORT',message='At least one RRT-by-SOFA cell empty');stop('Empty RRT/SOFA cell; author review required')}
contrasts<-data.table(contrast=c('MV_1_vs_0','RRT_1_vs_0','duration_48_vs_24_h','peak_0.30_vs_0.15'),variable=vars[c(1,2,6,5)],low=c(0,0,24,.15),high=c(1,1,48,.30))
support<-rbindlist(lapply(seq_len(nrow(contrasts)),function(j){v<-contrasts$variable[j];p<-quantile(dat[[v]],c(.05,.95),type=7);data.table(contrast=contrasts$contrast[j],variable=v,value=c(contrasts$low[j],contrasts$high[j]),p5=p[1],p95=p[2],within_p5_p95=c(contrasts$low[j],contrasts$high[j])>=p[1] & c(contrasts$low[j],contrasts$high[j])<=p[2])}))
fwrite(support,file.path(out,'contrast_support.csv'))
if(any(!support$within_p5_p95)) {status('STOP_SUPPORT',message='Preflight unexpected marginal support failure');stop('Unexpected support failure; author review required')}
contrast_x<-lapply(names(xs),function(m)lapply(seq_len(nrow(contrasts)),function(j){a<-copy(dat);b<-copy(dat);set(a,j=contrasts$variable[j],value=rep(contrasts$high[j],nrow(a)));set(b,j=contrasts$variable[j],value=rep(contrasts$low[j],nrow(b)));list(high=make_x(a,m=='sofa15'),low=make_x(b,m=='sofa15'))}));names(contrast_x)<-names(xs)
fit_model <- function(x,yy) {
  warns<-character();fit<-withCallingHandlers(glm.fit(x=x,y=yy,family=binomial(),control=glm.control(epsilon=1e-8,maxit=100)),warning=function(w){warns<<-c(warns,conditionMessage(w));invokeRestart('muffleWarning')})
  class(fit)<-c('glm','lm');fit$x<-x;fit$y<-yy
  bad<- !fit$converged || fit$rank!=ncol(x) || any(!is.finite(fit$coefficients))
  near<-any(grepl('numerically 0 or 1',warns)) || min(fit$fitted.values)<=1e-8 || max(fit$fitted.values)>=1-1e-8
  list(fit=fit,warnings=paste(warns,collapse=' | '),bad=bad,near=near)
}
summarize_contrasts <- function(beta,m) rbindlist(lapply(seq_len(nrow(contrasts)),function(j){cx<-contrast_x[[m]][[j]];ph<-plogis(drop(cx$high%*%beta));pl<-plogis(drop(cx$low%*%beta));ld<-drop((cx$high-cx$low)%*%beta);stopifnot(diff(range(ld))<1e-10);data.table(model=m,contrast=contrasts$contrast[j],standardized_risk_low=mean(pl),standardized_risk_high=mean(ph),risk_difference_pp=100*(mean(ph)-mean(pl)),adjusted_OR=exp(ld[1]))}))
fits<-list();main_coef<-list();main_con<-list();main_qc<-list();influence<-list()
for(m in names(xs)) {
 f<-fit_model(xs[[m]],y);fits[[m]]<-f$fit
 se<-sqrt(diag(summary(f$fit)$cov.unscaled));main_qc[[m]]<-data.table(model=m,n_events=length(y),n_patients=uniqueN(dat$subject_id),restarts=sum(y),predictor_df=ncol(xs[[m]])-1L,rank=f$fit$rank,converged=f$fit$converged,iterations=f$fit$iter,min_pred=min(f$fit$fitted.values),max_pred=max(f$fit$fitted.values),max_standard_error=max(se),warnings=f$warnings,near_separation_flag=f$near)
 if(f$bad||f$near) {fwrite(rbindlist(main_qc),file.path(out,'main_fit_checks.csv'));status('STOP_MAIN_FIT',model=m);stop('Main fit failure/near-separation warning')}
 main_coef[[m]]<-data.table(model=m,term=names(f$fit$coefficients),coefficient=as.numeric(f$fit$coefficients))
 main_con[[m]]<-summarize_contrasts(f$fit$coefficients,m)
 cd<-cooks.distance(f$fit);influence[[m]]<-data.table(model=m,rank=1:10,cooks_distance=sort(cd,decreasing=TRUE)[1:10])
}
fwrite(rbindlist(main_qc),file.path(out,'main_fit_checks.csv'))
fwrite(rbindlist(main_coef),file.path(out,'coefficients_point.csv'))
fwrite(rbindlist(main_con),file.path(out,'standardized_contrasts_point.csv'))
fwrite(rbindlist(influence),file.path(out,'cooks_top10_no_identifiers.csv'))
saveRDS(list(models=fits,knots=knots,contrasts=contrasts,standardization_target='original_1958_events',seed=seed),file.path(out,'restricted_work/main_models.rds'))
RNGkind('Mersenne-Twister','Inversion','Rejection');set.seed(seed)
groups<-split(seq_len(nrow(dat)),dat$subject_id)
bootcoef<-list();bootcon<-list();bootqc<-list();failures<-list();nf<-0L
for(b in 1:1000) {
 sampled<-sample.int(length(groups),length(groups),replace=TRUE);ix<-unlist(groups[sampled],use.names=FALSE)
 for(m in names(xs)) {
  key<-paste(b,m,sep='_')
  f<-tryCatch(fit_model(xs[[m]][ix,,drop=FALSE],y[ix]),error=function(e)list(error=conditionMessage(e)))
  if(!is.null(f$error)) {bad<-TRUE;near<-FALSE;why<-f$error;conv<-FALSE;iters<-NA_integer_;minp<-maxp<-NA_real_;rank<-NA_integer_} else {bad<-f$bad;near<-f$near;why<-f$warnings;conv<-f$fit$converged;iters<-f$fit$iter;minp<-min(f$fit$fitted.values);maxp<-max(f$fit$fitted.values);rank<-f$fit$rank}
  bootqc[[key]]<-data.table(replicate=b,model=m,converged=conv,rank=rank,iterations=iters,bootstrap_events=length(ix),restarts=sum(y[ix]),min_pred=minp,max_pred=maxp,near_separation_flag=near,success=!bad&&!near,warnings=why)
  if(bad||near) {nf<-nf+1L;failures[[key]]<-data.table(replicate=b,model=m,reason=if(near)paste('Near-separation:',why) else why);fwrite(rbindlist(failures),file.path(out,'bootstrap_failures.csv'));if(near||nf>10L){fwrite(rbindlist(bootqc),file.path(out,'bootstrap_fit_checks.csv'));status('STOP_BOOTSTRAP',replicate=b,failures=nf);stop('Bootstrap separation warning or failure threshold exceeded')}} else {
   bootcoef[[key]]<-data.table(replicate=b,model=m,term=names(f$fit$coefficients),coefficient=as.numeric(f$fit$coefficients))
   z<-summarize_contrasts(f$fit$coefficients,m);z[,replicate:=b];bootcon[[key]]<-z
  }
 }
 if(b%%25==0L) {logline('Bootstrap',b,'of1000; failed model fits',nf);fwrite(rbindlist(bootqc),file.path(out,'bootstrap_fit_checks.csv'));status('BOOTSTRAP_RUNNING',completed=b,failed_model_fits=nf)}
}
bc<-rbindlist(bootcoef);bt<-rbindlist(bootcon);bq<-rbindlist(bootqc)
fwrite(bc,file.path(out,'bootstrap_coefficients.csv'));fwrite(bt,file.path(out,'bootstrap_standardized_contrasts.csv'));fwrite(bq,file.path(out,'bootstrap_fit_checks.csv'))
fwrite(if(length(failures)) rbindlist(failures) else data.table(replicate=integer(),model=character(),reason=character()),file.path(out,'bootstrap_failures.csv'))
q<-function(x)quantile(x,c(.025,.975),type=7,names=FALSE)
cc<-bc[,.(coefficient_lower=q(coefficient)[1],coefficient_upper=q(coefficient)[2],successful_bootstrap=.N),by=.(model,term)]
cc<-merge(rbindlist(main_coef),cc,by=c('model','term'),sort=FALSE);fwrite(cc,file.path(out,'coefficients_with_bootstrap_ci.csv'))
long<-melt(bt,id.vars=c('model','contrast','replicate'),measure.vars=c('standardized_risk_low','standardized_risk_high','risk_difference_pp','adjusted_OR'),variable.name='measure')
cis<-long[,.(lower=q(value)[1],upper=q(value)[2],successful_bootstrap=.N),by=.(model,contrast,measure)]
pts<-melt(rbindlist(main_con),id.vars=c('model','contrast'),variable.name='measure',value.name='estimate')
fin<-merge(pts,cis,by=c('model','contrast','measure'),sort=FALSE);fwrite(fin,file.path(out,'standardized_contrasts_with_ci.csv'))
capture.output(sessionInfo(),file=file.path(out,'R_session_info.txt'))
meta<-list(seed=seed,repetitions=1000,failed_model_fits=nf,standardization_target_n=1958,elapsed_seconds=as.numeric(difftime(Sys.time(),start,units='secs')),quantile_type=7,interpretation='specified_adjustment_association_differences_not_causal',main_fit_count=2,bootstrap_fit_count=nrow(bq))
write_json(meta,file.path(out,'execution_summary.json'),auto_unbox=TRUE,pretty=TRUE,digits=NA)
status('COMPUTATION_COMPLETE',completed=1000,failed_model_fits=nf)
logline('Association computation finished')
