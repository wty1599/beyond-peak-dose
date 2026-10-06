# Repository disclosure-protected derivatives, drawn only from released summaries.
library(ggplot2)
library(patchwork)
source('code/00_setup/theme_manuscript.R')
read_summary <- function(name) data.table::fread(file.path('results/aggregate',name),encoding='UTF-8',data.table=FALSE)
db_levels <- names(cohort_colors)
percent_scale <- function(top=60) scale_y_continuous(limits=c(0,top),expand=expansion(mult=c(0,.02)))

# Structural flow only. Linked count families are withheld, including totals.
flow <- read_summary('Figure_1_structure.csv')
flow <- flow[!grepl('exclusions|exclusion',flow$stage,ignore.case=TRUE),]
flow$x <- match(flow$database,db_levels)
flow$y <- ave(seq_len(nrow(flow)),flow$database,FUN=function(z) 7-seq_along(z))
flow$label <- paste0(flow$stage,'\nCounts withheld')
flow$label <- gsub('First qualifying withdrawal','First qualifying\nwithdrawal',flow$label)
arrows <- flow[flow$y>1,]
p1 <- ggplot(flow,aes(x,y))+
  geom_segment(data=arrows,aes(xend=x,y=y-.42,yend=y-.59),linewidth=.25,
               arrow=arrow(length=unit(1.2,'mm')))+
  geom_rect(aes(xmin=x-.44,xmax=x+.44,ymin=y-.40,ymax=y+.40,color=database),fill='white',linewidth=.4)+
  geom_text(aes(label=label),size=2.7,lineheight=1.05)+
  scale_color_manual(values=cohort_colors,guide='none')+
  scale_x_continuous(breaks=1:3,labels=db_levels,limits=c(.45,3.55),position='top')+
  scale_y_continuous(limits=c(.4,6.7),expand=c(0,0))+
  labs(title='Cohort flow structure',subtitle='Public version: linked counts and exclusion counts withheld',
       caption='Steps follow the sequence used within each database; full flow counts are in the manuscript.')+
  theme_manuscript()+theme(axis.title=element_blank(),axis.line=element_blank(),axis.ticks=element_blank(),
                          axis.text.y=element_blank(),axis.text.x=element_text(face='bold',size=8))
write_pdf(p1,'Figure_1_public',155)

overall <- read_summary('Figure_2_fixed_times.csv')
point <- overall[overall$record_type%in%c('interval','point'),]
point$time <- factor(point$time_hours,levels=c(12,24,48,72))
pd <- position_dodge(width=.55)
p2a <- ggplot(point,aes(time,100*estimate,color=database,shape=database,group=database))+
  geom_errorbar(aes(ymin=100*lower,ymax=100*upper),position=pd,width=.13,linewidth=.4,na.rm=TRUE)+
  geom_point(position=pd,size=1.8)+scale_color_manual(values=cohort_colors)+scale_shape_manual(values=cohort_shapes)+
  percent_scale(45)+labs(title='A  Restart by the specified hour',x='Hours after cessation',y='Restart (%)')+theme_manuscript()
remaining <- overall[overall$record_type=='conditional_24_to_72h',]
remaining$database <- factor(remaining$database,levels=rev(db_levels))
p2b <- ggplot(remaining,aes(100*estimate,database,color=database,shape=database))+
  geom_errorbar(aes(xmin=100*lower,xmax=100*upper),orientation='y',width=.15,linewidth=.4)+
  geom_point(size=1.9)+scale_color_manual(values=cohort_colors,guide='none')+scale_shape_manual(values=cohort_shapes,guide='none')+
  scale_x_continuous(limits=c(0,20),breaks=seq(0,20,5))+labs(title='B  Remaining risk from 24 to 72 hours',
  subtitle='Alive, still in care and restart-free at 24 hours',x='Restart (%)',y=NULL)+theme_manuscript()
write_pdf((p2a/p2b)+plot_layout(heights=c(1.7,1))+
          plot_annotation(caption='Public fixed-time summaries; no interpolation. MIMIC-IV and SICdb: competing-event estimates.\nAmsterdamUMCdb: fixed-window proportions. Intervals are the available 95% bootstrap intervals.'),
          'Figure_2_public',170)

support <- read_summary('Figure_3_fixed_times.csv')
support <- support[support$record_type=='interval',]
support$time <- factor(support$time_hours,levels=c(12,24,72))
support$status <- ifelse(support$group=='MV=1','Support recorded','No support recorded')
p3 <- ggplot(support,aes(time,100*estimate,color=status,shape=status,group=status))+
  geom_errorbar(aes(ymin=100*lower,ymax=100*upper),position=pd,width=.12,linewidth=.4)+
  geom_point(position=pd,size=1.8)+facet_wrap(~database,nrow=1)+
  scale_color_manual(values=c('Support recorded'='#1B3A5C','No support recorded'='#777777'))+
  scale_shape_manual(values=c('Support recorded'=16,'No support recorded'=1))+
  percent_scale(60)+labs(title='Restart by ventilation status at cessation',
  subtitle='Public fixed-time summaries',x='Hours after cessation',y='Cumulative restart (%)',
  caption='Death and normal care-period ending are competing events. SICdb ventilation is an airway/device proxy.\nError bars: available 95% bootstrap intervals. No connecting lines or interpolated values are used.')+theme_manuscript()
write_pdf(p3,'Figure_3_public',100)

effects <- read_summary('Figure_4_source_data.csv')
stopifnot(nrow(effects)==9L,all(is.finite(effects$estimate)),sum(effects$panel=='A')==6L)
effects$label <- c('Mechanical ventilation\nMIMIC-IV','Ventilation\nSICdb (airway/device signal)',
                  'Ventilation\nAmsterdamUMCdb (recorded process)','Renal replacement therapy\nMIMIC-IV',
                  'Renal support\nSICdb (recorded CRRT, preceding 6 h)',
                  'Renal support\nAmsterdamUMCdb (recorded CVVH process)',db_levels)
effects$RD_display <- gsub('\u2212','-',effects$RD_display,fixed=TRUE,useBytes=TRUE)
forest <- function(d,title,xlim,breaks,axis_label,numbers=TRUE) {
  d$y <- rev(seq_len(nrow(d)))
  p <- ggplot(d,aes(estimate,y,color=database,shape=database))+
    geom_vline(xintercept=0,linewidth=.25,color='#777777',linetype=2)+
    geom_errorbar(aes(xmin=lower,xmax=upper),orientation='y',width=.18,linewidth=.45)+
    geom_point(size=1.8)+scale_color_manual(values=cohort_colors,guide='none')+
    scale_shape_manual(values=cohort_shapes,guide='none')+
    scale_x_continuous(limits=xlim,breaks=breaks)+scale_y_continuous(breaks=d$y,labels=d$label,
    limits=c(.4,max(d$y)+.6))+labs(title=title,x=axis_label,y=NULL)+theme_manuscript()
  textcol <- ggplot(d,aes(y=y))+geom_text(aes(x=0,label=RD_display),hjust=0,size=2.45)+
    {if(numbers) geom_text(aes(x=1.35,label=OR_display),hjust=0,size=2.45) else NULL}+
    scale_x_continuous(limits=c(0,if(numbers)2.75 else 1.45),expand=c(0,0))+
    scale_y_continuous(limits=c(.4,max(d$y)+.6))+
    labs(title=if(numbers)'RD, pp (95% CI)                 OR (95% CI)' else 'Change, pp (95% CI)')+
    theme_void(base_family='Arial')+theme(plot.title=element_text(size=7),plot.margin=margin(6,0,25,2))
  (p+textcol)+plot_layout(widths=if(numbers)c(1.4,1.25)else c(1.4,1.25))
}
p4a<-forest(effects[effects$panel=='A',],'A  Adjusted associations',c(-10,30),seq(-10,30,10),'Standardized risk difference (percentage points)')
p4b<-forest(effects[effects$panel=='B',],'B  Change after adding dose category',c(-2.5,2.5),seq(-2,2,1),'Change in ventilation risk difference (percentage points)',FALSE)
write_pdf((p4a/p4b)+plot_layout(heights=c(1.7,1))+
          plot_annotation(caption='Panel A uses each cohort-specific primary adjusted model. Panel B uses the common model form\nwith and without peak-dose category on the same sample within each cohort. Estimates are not pooled.'),
          'Figure_4_public',180)

dose <- read_summary('Figure_S6_fixed_times.csv')
dose$time <- factor(dose$time_hours,levels=c(12,24,48,72))
p6 <- ggplot(dose,aes(time,100*estimate,color=group,shape=group,group=group))+
  geom_errorbar(aes(ymin=100*lower95,ymax=100*upper95),position=pd,width=.12,linewidth=.4,na.rm=TRUE)+
  geom_point(position=pd,size=1.8)+facet_wrap(~database,nrow=1)+
  scale_color_manual(values=c('<0.20'='#1B3A5C','0.20-<0.50'='#527BA4','>=0.50'='#C1662F'))+
  scale_shape_manual(values=c('<0.20'=15,'0.20-<0.50'=16,'>=0.50'=17))+
  percent_scale(60)+labs(title='Restart by peak norepinephrine-equivalent dose',
  subtitle='Dose in micrograms/kg/min; public fixed-time summaries',x='Hours after cessation',y='Cumulative restart (%)',
  caption='Death and normal care-period ending are competing events. Available 95% bootstrap intervals are shown.\nNo connecting lines or interpolated values are used.')+theme_manuscript()
write_pdf(p6,'Figure_S6_public',100)
