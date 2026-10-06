suppressPackageStartupMessages(library(data.table))
# Weighted AJ implementation mechanically retained from the approved prior analysis.
# All participants enter at q=4h; no estimated curves before q are delivered.
aj <- function(time,status,weight=rep(1,length(time))) {
 stopifnot(length(time)>0,length(time)==length(status),length(time)==length(weight),all(is.finite(time)),all(weight>=0))
 u<-sort(unique(time));j<-match(time,u)
 tab<-function(v){z<-numeric(length(u));z0<-rowsum(v,j,reorder=FALSE);z[as.integer(rownames(z0))]<-z0[,1];z}
 nt<-tab(weight);nr<-rev(cumsum(rev(nt)));dd<-tab(weight*(status>0));d1<-tab(weight*(status==1));d2<-tab(weight*(status==2));d3<-tab(weight*(status==3))
 hazard<-ifelse(nr>0,dd/nr,0);safter<-cumprod(1-hazard);sbefore<-c(1,head(safter,-1L))
 data.table(time_hours=u,restart=cumsum(sbefore*ifelse(nr>0,d1/nr,0)),death=cumsum(sbefore*ifelse(nr>0,d2/nr,0)),
  care_end=cumsum(sbefore*ifelse(nr>0,d3/nr,0)),event_free=safter,at_risk_before=nr)
}
stepat<-function(d,t,col){j<-findInterval(t,d$time_hours);c(if(col=='event_free')1 else 0,d[[col]])[j+1L]}
# Independent event-by-event loop matching the original MIMIC sap5_7/S20 math.
aj_reference<-function(time,status,grid,weight=rep(1,length(time)),cause=1L){
 vapply(grid,function(h){surv<-1;cif<-0
  for(t in sort(unique(time[status>0&time<=h&weight>0]))){nr<-sum(weight[time>=t]);dd<-sum(weight[time==t&status>0]);dc<-sum(weight[time==t&status==cause]);if(nr>0){cif<-cif+surv*dc/nr;surv<-surv*(1-dd/nr)}}
  cif},numeric(1))
}
