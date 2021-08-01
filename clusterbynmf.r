#################################################
#  File Name:cluster.r
#  Author: xingpengwei
#  Mail: xingwei421@qq.com
#  Created Time: Thu 07 Jan 2021 01:23:41 PM UTC
#################################################

library(NMF)
library(foreach)
library(doParallel)
library(bigmemory)
library(synchronicity)
library(limma)
library(edgeR)
args <- commandArgs(TRUE)
all_data = read.table(args[1],header=T)
gene_length=all_data[,3]-all_data[,2]
all_data=all_data[,4:length(all_data[1,])]
all_data2=apply(all_data,2,as.numeric)
data_norm=rpkm(all_data2,gene.length = gene_length)
mysum=apply(data_norm,1,sum)
mysum = as.data.frame(mysum)
data_norm2=data_norm[which(mysum>0),]
res2 <- nmf(data_norm2, 2, 'nsNMF', seed=123456,nrun=100,.pbackend='par',.opt='vp50')
res3 <- nmf(data_norm2, 3, 'nsNMF', seed=123456,nrun=100,.pbackend='par',.opt='vp50')
res4 <- nmf(data_norm2, 4, 'nsNMF', seed=123456,nrun=100,.pbackend='par',.opt='vp50')
res5 <- nmf(data_norm2, 5, 'nsNMF', seed=123456,nrun=100,.pbackend='par',.opt='vp50')
res6 <- nmf(data_norm2, 6, 'nsNMF', seed=123456,nrun=100,.pbackend='par',.opt='vp50')
res7 <- nmf(data_norm2, 7, 'nsNMF', seed=123456,nrun=100,.pbackend='par',.opt='vp50')
res8 <- nmf(data_norm2, 8, 'nsNMF', seed=123456,nrun=100,.pbackend='par',.opt='vp50')
pdf("cluster.rank2.nmf.pdf", width=10, height=10)
consensusmap(res2)
dev.off()
pdf("cluster.rank3.nmf.pdf", width=10, height=10)
consensusmap(res3)
dev.off()
pdf("cluster.rank4.nmf.pdf", width=10, height=10)
consensusmap(res4)
dev.off()
pdf("cluster.rank5.nmf.pdf", width=10, height=10)
consensusmap(res5)
dev.off()
pdf("cluster.rank6.nmf.pdf", width=10, height=10)
consensusmap(res6)
dev.off()
pdf("cluster.rank7.nmf.pdf", width=10, height=10)
consensusmap(res7)
dev.off()
pdf("cluster.rank8.nmf.pdf", width=10, height=10)
consensusmap(res8)
dev.off()
save.image("my.Rdata")
