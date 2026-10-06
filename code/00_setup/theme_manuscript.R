library(ggplot2)
cohort_colors <- c('MIMIC-IV'='#1B3A5C','SICdb'='#C1662F','AmsterdamUMCdb'='#666666')
cohort_shapes <- c('MIMIC-IV'=15,'SICdb'=16,'AmsterdamUMCdb'=17)
theme_manuscript <- function() {
  theme_classic(base_size=8,base_family='Arial') +
    theme(axis.text=element_text(size=7,color='black'),
          axis.title=element_text(size=8),
          axis.line=element_line(linewidth=.18),
          axis.ticks=element_line(linewidth=.18),
          legend.position='bottom',legend.text=element_text(size=7),
          legend.title=element_blank(),
          strip.background=element_blank(),strip.text=element_text(size=8,face='bold'),
          plot.title=element_text(size=9,face='bold'),
          plot.subtitle=element_text(size=7),
          plot.caption=element_text(size=7,hjust=0),
          plot.margin=margin(6,9,6,6))
}
write_pdf <- function(plot,name,height) {
  dir.create('figures',showWarnings=FALSE)
  ggsave(file.path('figures',paste0(name,'.pdf')),plot,width=180,height=height,
         units='mm',device=cairo_pdf,bg='white')
}
