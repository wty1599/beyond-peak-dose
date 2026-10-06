# Convex feasibility screen, not an alternative association estimator.
separation_check <- function(x,y,model){
  scales<-apply(abs(x),2,max);scales[scales==0]<-1
  z<-sweep(x,2,scales,'/')*(2*y-1);p<-ncol(z)
  obj<-c(colSums(z),-colSums(z))
  ans<-boot::simplex(obj,A1=rbind(-cbind(z,-z),rep(1,2*p)),b1=c(rep(0,nrow(z)),1),maxi=TRUE,n.iter=50000,eps=1e-10)
  beta<-ans$soln[seq_len(p)]-ans$soln[p+seq_len(p)]
  margin<-drop(z%*%beta);value<-sum(margin)
  sep<-is.finite(value)&&value>1e-7&&min(margin)>-1e-7
  data.table(model=model,LP_solved=ans$solved,objective=value,minimum_margin=min(margin),maximum_margin=max(margin),separation_or_quasi_separation=sep,pass=ans$solved==1L&&!sep)
}
tests<-list(overlap=list(x=c(0,1,0,1),y=c(0,0,1,1),expected=FALSE),complete=list(x=c(-2,-1,1,2),y=c(0,0,1,1),expected=TRUE),quasi=list(x=c(-1,0,0,1),y=c(0,0,1,1),expected=TRUE))
synthetic<-rbindlist(lapply(names(tests),function(nm){q<-tests[[nm]];z<-separation_check(cbind(1,q$x),q$y,nm);z[,expected:=q$expected];z}))
stopifnot(all(synthetic$LP_solved==1L),all(synthetic$separation_or_quasi_separation==synthetic$expected))
fwrite(synthetic,file.path(a,'separation_synthetic_unit_tests.csv'))
