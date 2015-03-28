#!/bin/env Rscript

library(dplyr)
library(ggplot2)

gender <- rep(c('Men','Women'),5)
age <- c(rep('15-49',2),rep('50-59',2),rep('60-69',2),
	rep('70-79',2),rep('80-99',2))
survival <- c(17,23,
		9,10,
		6,6,
		4,4,
		2,3)
p <- data.frame(gender,age,survival) %>%
	ggplot(aes(x=age,fill=gender,y=survival))+
		geom_bar(stat='identity',position='dodge')+
		labs(x='Age group',y='5-year survival rate (%)',
			fill= ' ')+
		theme(text = element_text(color='black'))
ggsave(p,file='survivalRate.pdf')
