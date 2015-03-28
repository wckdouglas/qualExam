#!/usr/bin/env Rscript

library(dplyr)
library(tidyr)
library(ggplot2)

a <- c( c(55,32,13),c(30,32,13))
b <- rep(c('CTC EMT+','CTC platelet','CTC epithelial'),2)
d <- c(rep('Published',3),rep('Expected',3))
p <- data.frame(frac = a,name = b,annotate=d) %>%
	spread(annotate,frac) %>%
	mutate(Expected = Expected/sum(Expected)*100,
		Published = Published/sum(Published)*100) %>%
	gather(annotate,frac,-name) %>%
	mutate(annotate = factor(annotate,levels=c('Published','Expected')),
		name = factor(name,levels=(c('CTC EMT+','CTC platelet','CTC epithelial')))) %>%
	filter(annotate=='Published') %>%
	rev() %>%
	ggplot(aes(x='CTC Heterogenity',y=frac,fill=name,color=name,
			order = -as.numeric(name)))+
		geom_bar(stat='identity',alpha=0.8)+
#		facet_wrap(~annotate)+
		labs(x= ' ',y='Count',fill='Subtype',color='Subtype')+
		theme(panel.background=element_blank(),
			panel.grid=element_blank(),
			axis.text.x=element_blank())
ggsave(p,file='CTC_heterogeneity.pdf',height=7,width=4)
