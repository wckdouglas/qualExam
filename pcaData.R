#!/usr/bin/env Rscript

library(data.table)
library(dplyr)
library(reshape2)
library(ggplot2)
library(tidyr)
library(grid)

annotate <- function(name){
	if (grepl('^TuGMP',name)){
		return('primary tumor (TuGMP3)')
	}else if (grepl('MEF',name)){
		return('mouse embryonic fibroblasts')
	}else if (grepl('^MP',name)){
		return('circulating tumor cells')
	}else if (grepl('WBC',name)){
		return('leukocytes WBCs')
	}else if (grepl('nb',name)){
		return('mouse PDAC cell line (NB508)')}
}

setwd('~/qualexam')
genes <- c('Met','Snai2','Bcl2','Vegfa','Mmp14','Ldha',
		'Pdk1','Twist1','Igf2','Il10','Mif','Vim')
dat <- fread('mouseCTC_DESeq.tsv') %>%
		setnames(make.names(names(.))) %>%
		filter(symbol %in% genes) %>%
		select(4,5,7:193) %>%
		gather(sample,counts,-symbol,-name) %>%
		filter(grep('^TuGMP|^MP',sample)) %>%
		select(-name) %>%
		spread(symbol,counts) %>%
		gather(genes,counts,-sample) %>%
		mutate(counts = ifelse(is.na(counts),0,counts)) %>%
		spread(genes,counts) 

pca <- prcomp(dat[,-1,with=F],scale=T) 

mat <- data.frame(pca$x,dat[,1,with=F])  %>%
	mutate(annot = sapply(sample,annotate))

rotation_data <- data.table(pca$rotation, variable=row.names(pca$rotation))

arrowColor = 'black'

arrow_style <- arrow(length = unit(0.2, "cm"),
                     type = "closed")
p <- ggplot(mat, aes(x=PC1,y=PC2))+
	geom_text(aes(label=sample,color=annot)) +
	geom_segment(data=rotation_data,
              aes(xend=PC1*4, yend=PC2*4),
              x=0, y=0,
	      arrow=arrow_style,col=arrowColor)+
	geom_text(data=rotation_data,
		aes(x=PC1*4,y=PC2*4,label=variable),
		hjust=0,col=arrowColor) +
	coord_fixed() + 
	labs (color='Cell type')

ggsave(p,file='pca_hypoxia_signaling.pdf',height=10,width=10)
