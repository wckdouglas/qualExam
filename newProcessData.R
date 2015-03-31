#!/usr/bin/env Rscript

library(data.table)
library(dplyr)
library(reshape2)
library(ggplot2)
library(tidyr)
setwd('~/qualexam')
annotationFile <- fread('sampleAnnotation.csv',header=F) %>%
			setnames(c('sample','annotation')) %>%
			mutate(sample = make.names(sample))
genes <- c('Hif1a','Met','Snai2','Bcl2','Vegfa','Mmp14','Ldha',
		'Pdk1','Twist1','Igf2','Il10','Mif','Vim')
dat <- fread('mouseCTC_DESeq.tsv') %>%
		filter(symbol %in% genes) %>%
		select(4,5,7:193) %>%
		gather(sample,counts,-symbol,-name) %>%
		inner_join(annotationFile) %>%
		filter(annotation !='primary tumor') %>%
		mutate(annotation = ifelse(annotation=='primary tumor 2',
						'primary tumor',annotation))

hif <-	dat %>% 
		select (-name,-sample) %>%
		group_by (symbol,annotation) %>%
		summarize(count=mean(counts),
			stdv=(sd(counts)),
			sampleSize=n()) %>%#/sqrt(sum(record))) 
		mutate(annotation = paste0(annotation,' (n=',sampleSize,')'))

Ttest <- dat %>%
		select(-sample) %>% 
		group_by(symbol) %>%
		summarize(p = t.test(.[annotation=='CTC-EMT']$counts,
					.[annotation=='primary tumor']$counts)$p.value,
			ymax = max(mean(counts)))



p <- ggplot(data=hif,aes(x=annotation,y=count))+
	geom_bar(aes(color=annotation,
			fill=annotation),
		alpha=0.8,
		stat='identity',
		width=0.9,
		position=position_dodge(width=0.7))+
#	geom_errorbar(data=hif,
#		aes(ymax=count+stdv,
#			ymin=count-stdv),
#		position=position_dodge(width=0.9),
#		width=0.25)+
	geom_text(data=Ttest,x=1.5,size=3,
		aes(y=ymax,label=paste('p-value =',signif(p,3))))+
	facet_wrap(~symbol,nrow=2,scales='free_y')+
	labs(x=' ',y='normalized mean count')+
	theme(axis.ticks.x=element_blank(),
		axis.text.x=element_blank(),
		strip.text=element_text(size=15),
		axis.title.y=element_text(size=15),
		legend.title=element_text(size=15),
		legend.text=element_text(size=14),
		legend.position = 'right')
ggsave(p,file='downstream_genes.pdf',width=14,height=8)
#ggsave(p,file='qualexam_HIF.pdf',width=7,height=7)
