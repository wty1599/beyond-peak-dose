options(stringsAsFactors=FALSE)
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(digest);library(jsonlite)})
source('code/00_setup/paths.R')
out<-private_output('mice_diagnostics')
source_file<-file.path(private_output('risk_model'),'sap3_main_mids_v1_3.rds')
source_hash <- digest(file=source_file,algo='sha256')
m <- readRDS(source_file)
stopifnot(inherits(m,'mids'),m$m==20L,m$iteration==10L)
targets <- names(m$method)[m$method!='']
stopifnot(length(targets)==4L, identical(dim(m$chainMean),c(14L,10L,20L)),identical(dim(m$chainVar),dim(m$chainMean)))
labels <- c(predictor_mean_map_4h='Mean MAP',predictor_map_source_invasive='MAP source (factor codes)',predictor_mean_heart_rate_4h='Mean heart rate',predictor_urine_output_24h_ml='Urine output')
trace <- rbindlist(lapply(targets,function(v) rbindlist(lapply(1:20,function(j) data.table(variable=v,chain=j,iteration=1:10,chain_mean=as.numeric(m$chainMean[v,,j]),chain_variance=as.numeric(m$chainVar[v,,j]))))))
stopifnot(nrow(trace)==800L,all(is.finite(trace$chain_mean)),all(is.finite(trace$chain_variance)))
fwrite(trace,file.path(out,'saved_chain_mean_variance.csv'))
long <- melt(trace,id.vars=c('variable','chain','iteration'),variable.name='statistic',value.name='value')
long[,variable_label:=factor(labels[variable],levels=labels[targets])]
long[,statistic:=factor(statistic,levels=c('chain_mean','chain_variance'),labels=c('Mean of imputed values','Variance of imputed values'))]
# Stable per-chain encoding: five colours x four line styles, never an outcome encoding.
palette <- rep(c('#1B3A5C','#C1662F','#6E9BC5','#7A7A7A','#333333'),4)
styles <- rep(c('solid','longdash','dotted','dotdash'),each=5)
theme_trace <- theme_classic(base_size=9,base_family='Arial')+theme(legend.position='bottom',legend.title=element_text(size=8),legend.text=element_text(size=7),strip.background=element_blank(),strip.text=element_text(size=9),axis.line=element_line(linewidth=.15),axis.ticks=element_line(linewidth=.15),plot.margin=margin(5,8,5,5))
p <- ggplot(long,aes(iteration,value,group=chain,colour=factor(chain),linetype=factor(chain)))+geom_line(linewidth=.3,alpha=.72)+
  facet_wrap(vars(variable_label,statistic),ncol=2,scales='free_y')+
  scale_x_continuous(breaks=c(1,2,4,6,8,10),limits=c(1,10))+
  scale_colour_manual(values=palette,name='Imputation chain')+scale_linetype_manual(values=styles,name='Imputation chain')+
  guides(colour=guide_legend(nrow=2,byrow=TRUE),linetype=guide_legend(nrow=2,byrow=TRUE))+
  labs(x='Iteration',y=NULL)+theme_trace
ggsave(file.path(out,'MICE_saved_traces_overview.pdf'),p,width=210,height=250,units='mm',device=cairo_pdf)
ggsave(file.path(out,'MICE_saved_traces_overview.png'),p,width=210,height=250,units='mm',dpi=180)
# Individual panels retain all chains for close inspection.
for(v in targets){
  pv <- ggplot(long[variable==v],aes(iteration,value,group=chain,colour=factor(chain),linetype=factor(chain)))+geom_line(linewidth=.4,alpha=.8)+
    facet_wrap(vars(statistic),ncol=2,scales='free_y')+scale_x_continuous(breaks=1:10,limits=c(1,10))+
    scale_colour_manual(values=palette,name='Imputation chain')+scale_linetype_manual(values=styles,name='Imputation chain')+
    guides(colour=guide_legend(nrow=2,byrow=TRUE),linetype=guide_legend(nrow=2,byrow=TRUE))+labs(x='Iteration',y=NULL,subtitle=labels[v])+theme_trace
  ggsave(file.path(out,paste0(v,'.png')),pv,width=210,height=100,units='mm',dpi=180)
}
stopifnot(identical(source_hash,digest(file=source_file,algo='sha256')))
write_json(list(source=source_file,sha256=source_hash,m=20,iterations=10,targets=targets,trace_rows=nrow(trace),imputation_reruns=0,all_trace_values_finite=TRUE,source_unchanged=TRUE,visual_convergence_judgement='PENDING_VISUAL_REVIEW'),file.path(out,'trace_extraction_qc.json'),pretty=TRUE,auto_unbox=TRUE)
cat('SAVED_TRACES_EXTRACTED_NO_IMPUTATION_RERUN\n')
