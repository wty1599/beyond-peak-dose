options(stringsAsFactors=FALSE,digits=17,warn=1)
suppressPackageStartupMessages({library(data.table);library(digest);library(jsonlite)})
source('code/00_setup/paths.R')
out<-private_output('mimic_burden');w<-file.path(out,'restricted_work')
start<-Sys.time()
paths<-c(cohort=release_path('mimic_selected_cohort'),intervals=file.path(w,'nee_after_restart.csv'))
before <- vapply(paths,function(p)digest(file=p,algo='sha256'),character(1))
parse_time <- function(x)as.numeric(as.POSIXct(x,tz='UTC'))
d <- fread(paths['cohort'])[restart_72h==1]
n <- fread(paths['intervals']);stopifnot(nrow(d)==680,uniqueN(d$stay_id)==680,setequal(d$stay_id,n$stay_id),all(is.finite(n$nee)),all(n$nee>0))
n[,a:=parse_time(raw_starttime)];n[,b:=parse_time(raw_endtime)]
stopifnot(all(n$b>n$a))
result <- rbindlist(lapply(seq_len(nrow(d)),function(i){
 r <- d[i];rr <- parse_time(r$cif_time);end <- parse_time(r$ascertained_followup_end);tt <- parse_time(r$t0)
 z <- n[stay_id==r$stay_id][order(a,b)]
 stopifnot(length(rr)==1,is.finite(end),end>rr,min(z$a)==rr,rr>=tt+4*3600,end<=tt+72*3600)
 z[,aa:=pmax(a,rr)];z[,bb:=pmin(b,end)];stopifnot(all(z$bb>z$aa))
 prior_end <- c(-Inf,head(cummax(z$bb),-1))
 hours <- sum(pmax(0,z$bb-pmax(z$aa,prior_end))/3600)
 cls <- if(max(z$b)>end)'POSITIVE_RECORD_CROSSES_OBSERVATION_END' else if(max(z$b)==end)'TERMINAL_TIE_COMPLETENESS_UNRESOLVED' else 'CESSATION_OBSERVED_BEFORE_WINDOW_END'
 stopifnot(hours>0,hours<=(end-rr)/3600+1e-10)
 data.table(stay_id=r$stay_id,subject_id=r$subject_id,first_restart=r$cif_time,followup_end=r$ascertained_followup_end,
   end_reason=r$followup_end_reason,completeness=cls,positive_hours=hours,peak_nee=max(z$nee),
   observable_hours=(end-rr)/3600,hours_to_restart=(rr-tt)/3600,positive_segments=nrow(z),
   overlap_segments=sum(z$aa<prior_end),source_ends_after_window=sum(z$b>end))
}))
fwrite(result,file.path(w,'restart_burden_event_records.csv'))
summary <- rbindlist(lapply(c('ALL',sort(unique(result$completeness))),function(g){
 z <- if(g=='ALL')result else result[completeness==g]
 rbindlist(lapply(c('positive_hours','peak_nee','observable_hours'),function(v){q <- quantile(z[[v]],c(.25,.5,.75),names=FALSE,type=7)
 data.table(group=g,measure=v,n=nrow(z),valid_n=sum(is.finite(z[[v]])),q1=q[1],median=q[2],q3=q[3],
 role=if(v=='observable_hours')'Observation-window metadata, not a third burden outcome' else 'Prespecified-in-addendum continuous descriptive burden')
 }))
}))
fwrite(summary,file.path(out,'burden_summary.csv'))
counts <- result[,.N,by=.(end_reason,completeness)][order(end_reason,completeness)]
fwrite(counts,file.path(out,'observation_end_and_completeness.csv'))
qc <- list(status='COMPUTED_PENDING_INDEPENDENT_REVIEW',existing_restart_n=nrow(d),burden_n=nrow(result),missing_burden=0,
  unidentified_followup_end=0,first_positive_start_mismatches=0,overlapping_source_segments=sum(result$overlap_segments),
  source_intervals=nrow(n),model_fits=0,tests=0,outcome_changes=0,late_restarts_not_assumed_lower_burden=TRUE,
  no_new_categorical_duration_threshold=TRUE,terminal_ties_not_complete=TRUE,elapsed_seconds=as.numeric(difftime(Sys.time(),start,units='secs')))
write_json(qc,file.path(out,'BURDEN_QC.json'),pretty=TRUE,auto_unbox=TRUE)
after <- vapply(paths,function(p)digest(file=p,algo='sha256'),character(1));stopifnot(identical(before,after))
fwrite(data.table(role=names(paths),path=unname(paths),before=unname(before),after=unname(after),unchanged=before==after),file.path(out,'source_protection.csv'))
cat('Seed not applicable; deterministic descriptive calculation; source query was READ ONLY and rolled back.\n')
print(summary);print(counts);print(qc);print(sessionInfo());
