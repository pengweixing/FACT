#################################################
#  File Name:process.r
#  Author: xingpengwei
#  Mail: xingwei421@qq.com
#  Created Time: Fri 25 Jun 2021 03:30:09 PM UTC
#################################################

library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(Organism.dplyr)

data = read.table("gene_expression_matrix",header=T)
###gene    RJ051C  RJ053C  RJ043C 
#  OR4F5   1	1	2
#  OR4F29  0	1	2
data2 = data
data3 = data2[,2:ncol(data2)]
sample = colnames(data3)
gene = data2[,1]
src <- src_ucsc("Homo sapiens")
tx2=as.data.frame(transcripts(src,columns = c("symbol")))
for(i in 1:ncol(data3)){
	temp = data3[,i]
	temp_q = quantile(temp)
	temp_q2=as.data.frame(t(temp_q))
	quan1_index = which(temp<temp_q[2])
	quan2_index = which(temp>=temp_q[2] & temp<temp_q[3])
	quan3_index = which(temp>=temp_q[3] & temp<temp_q[4])
	quan4_index = which(temp>=temp_q[4])
	quan1_gene = gene[quan1_index]
	quan2_gene = gene[quan2_index]
	quan3_gene = gene[quan3_index]
	quan4_gene = gene[quan4_index]
	index1 = which(!is.na(match(tx2$symbol,quan1_gene)))
	index2 = which(!is.na(match(tx2$symbol,quan2_gene)))
	index3 = which(!is.na(match(tx2$symbol,quan3_gene)))
	index4 = which(!is.na(match(tx2$symbol,quan4_gene)))
	quantile1_tx = tx2[index1,]
	quantile2_tx = tx2[index2,]
	quantile3_tx = tx2[index3,]
	quantile4_tx = tx2[index4,]
	quantile1_tx2 = quantile1_tx[,c(1,2,3,7,8,5)]
	quantile2_tx2 = quantile2_tx[,c(1,2,3,7,8,5)]
	quantile3_tx2 = quantile3_tx[,c(1,2,3,7,8,5)]
	quantile4_tx2 = quantile4_tx[,c(1,2,3,7,8,5)]
	if(!dir.exists(sample[i])){
	dir.create(sample[i])
	}
	write.table(quantile1_tx2,file=paste0(sample[i],"/quantile1_tx2.bed"),sep="\t",row.names=F,col.names=F,quote=F)
	write.table(quantile2_tx2,file=paste0(sample[i],"/quantile2_tx2.bed"),sep="\t",row.names=F,col.names=F,quote=F)
	write.table(quantile3_tx2,file=paste0(sample[i],"/quantile3_tx2.bed"),sep="\t",row.names=F,col.names=F,quote=F)
	write.table(quantile4_tx2,file=paste0(sample[i],"/quantile4_tx2.bed"),sep="\t",row.names=F,col.names=F,quote=F)

	write.table(quan1_gene,file=paste0(sample[i],"/quantile1_gene.bed"),sep="\t",row.names=F,col.names=F,quote=F)
	write.table(quan2_gene,file=paste0(sample[i],"/quantile2_gene.bed"),sep="\t",row.names=F,col.names=F,quote=F)
	write.table(quan3_gene,file=paste0(sample[i],"/quantile3_gene.bed"),sep="\t",row.names=F,col.names=F,quote=F)
	write.table(quan4_gene,file=paste0(sample[i],"/quantile4_gene.bed"),sep="\t",row.names=F,col.names=F,quote=F)
}
