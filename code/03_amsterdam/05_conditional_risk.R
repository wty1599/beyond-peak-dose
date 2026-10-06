options(stringsAsFactors=FALSE,warn=1,digits=17)
Sys.setenv(OMP_NUM_THREADS=1,OPENBLAS_NUM_THREADS=1,MKL_NUM_THREADS=1)
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(Hmisc))
suppressPackageStartupMessages(library(jsonlite))
setDTthreads(1)
source('code/00_setup/paths.R')
stage<-private_output('amsterdam_descriptive')
started<-Sys.time()

# Requested A1 fixed-window proportion, no CIF and no new endpoint extraction.
a<-fread(file.path(stage,'a1_rows.csv'))
target<-a[condition=='eligible']; valid<-target[restart72 %in% c('YES','NO')]
groups<-split(seq_len(nrow(target)),target$patientid)
set.seed(202610021L); B<-1000L
draws<-numeric(B)
for(b in seq_len(B)) {
  ix<-unlist(groups[sample.int(length(groups),length(groups),replace=TRUE)],use.names=FALSE)
  z<-target[ix];z<-z[restart72 %in% c('YES','NO')]
  draws[b]<-100*mean(z$restart72=='YES')
}
lim<-quantile(draws,c(.025,.975),type=7,names=FALSE)
fwrite(data.table(estimate=100*mean(valid$restart72=='YES'),lower=lim[1],upper=lim[2],
                 condition_n=nrow(target),n=sum(valid$restart72=='YES'),N=nrow(valid),
                 seed=202610021L,bootstrap_repetitions=B),file.path(stage,'a1_statistic.csv'))
hours<-fread(file.path(stage,'a4_duration.csv'))$hours
q<-quantile(hours,c(.25,.5,.75),type=7,names=FALSE,na.rm=TRUE)
fwrite(data.table(N=sum(!is.na(hours)),estimate=q[2],lower=q[1],upper=q[3]),file.path(stage,'a4_duration_summary.csv'))
cat('A1 fixed-window proportion and A4 duration completed.\n')

