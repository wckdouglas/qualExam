#!/usr/bin/env Rscript

#require(data.table)
#suppressMessages(require(DESeq2))
#require(BiocParallel)
#register(MulticoreParam(12))

setwd('~/qualexam')
#dat <- fread('mouseCTC.tsv')
#mat <- dat[,7:193,with=F]
#setnames(mat,
#	make.names(colnames(mat),
#		unique=T))

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
#
#cells <- sapply(colnames(mat),annotate)
#
#colData <- data.table(names=colnames(mat),
#			cellType=factor(cells))
#dds <- DESeqDataSetFromMatrix(countData=mat,
#			colData = colData,design= ~cellType)
#
#dds <- DESeq(dds,
#		betaPrior=FALSE)
#
#dat <- data.table(dat[,1:6,with=F],counts(dds,normalized=T))
#setnames(dat,c(colnames(dat[,1:6,with=F]),colnames(mat)))
#write.table(dat,'mouseCTC_DESeq.tsv',sep='\t',quote=F,row.names=F,col.names=T)

#======================================

require(data.table)
require(dplyr)
require(reshape2)
require(ggplot2)
require(tidyr)
setwd('~/qualexam')
dat <- fread('mouseCTC_DESeq.tsv') %>%
		filter(grep('Hif1|Vegf|Snai|Zeb|Arnt',symbol)) %>%
		select(4,5,7:193) %>%
		gather(sample,counts,-symbol,-name) %>%
		mutate(record = 1,
			annotation=sapply(sample,annotate)) 
hif <-	dat %>% 
		select (c(symbol,annotation,counts,record)) %>%
		filter (annotation %in% c('primary tumor (TuGMP3)','circulating tumor cells' )) %>%
		group_by (annotation,symbol) %>%
		summarize(count=mean(counts),stdv=(sd(counts)))#/sqrt(sum(record)))
Ttest <- dat %>%
		spread(annotation,counts) %>%
		select(1:5,9) %>%
		setnames(make.names(names(.))) %>%
		mutate(circulating.tumor.cells = ifelse(is.na(circulating.tumor.cells),0,
								circulating.tumor.cells),
			primary.tumor..TuGMP3. = ifelse(is.na(primary.tumor..TuGMP3.),0,
								primary.tumor..TuGMP3.)) %>%
		group_by(symbol) %>%
		summarize(pvalue=t.test(as.numeric(circulating.tumor.cells),
					as.numeric(primary.tumor..TuGMP3.))$p.value)


#t.test(hif[symbol=='Hif1a'][annotation=='primary tumor (TuGMP3)']$counts,hif[symbol=='Hif1a'][annotation=='primary tumor (TuGMP3)']$counts)

p <- ggplot(data=hif,aes(x=annotation,y=count))+
	geom_bar(aes(color=annotation,
			fill=annotation),
		stat='identity',
		width=0.9,
		position=position_dodge(width=0.9))+
#	geom_errorbar(data=hif,
#		aes(ymax=count+stdv,
#			ymin=count-stdv),
#		position=position_dodge(width=0.9),
#		width=0.25)+
	geom_text(data=Ttest,x=1.5,y=450,
		aes(label=paste('p-value =',signif(pvalue,3))))+
	facet_wrap(~symbol,ncol=5)+
	labs(x=' ',y='normalized mean count')+
	theme(axis.text.x=element_blank(),
		axis.ticks.x=element_blank(),
		strip.text=element_text(size=15))
ggsave(p,file='qualexam_HIF.pdf',width=10,height=5)
