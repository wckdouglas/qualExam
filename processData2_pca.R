#!/bin/env Rscript

require(data.table)
require(dplyr)
require(tidyr)
require(ggplot2)

#====== set up sample  annotation=======
annotate <- function(name){
	if (grepl('GMP',name)){
		return('primary tumor (TuGMP3)')
	}else if (grepl('MEF',name)){
		return('mouse embryonic fibroblasts')
	}else if (grepl('MP',name)){
		return('circulating tumor cells')
	}else if (grepl('WBC',name)){
		return('leukocytes WBCs')
	}else if (grepl('nb',name)){
		return('mouse PDAC cell line (NB508)')}
}

#=========read data 
setwd('~/qualExam')
dat <- fread('mouseCTC_DESeq.tsv') %>% 
	filter(grepl('Hif1|Snai|Vegf',symbol)) %>% 
	select(-1:-3,-5:-6) %>% 
	gather(sample,count,-symbol) %>%
	spread(symbol,count) 


# run pca
pca <- dat %>%
	select(-sample) %>%
	prcomp(center = TRUE, scale. = TRUE)

plotTable <- data.table(sample=dat$sample,pca$x) %>%
		mutate(annot = sapply(sample,annotate),
			sampleName = make.names(annot))

ggplot(data=plotTable,aes(x=PC1,y=PC2,label=sampleName)) +
	geom_text(aes(color=annot))+ coord_fixed()


