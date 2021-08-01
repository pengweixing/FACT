#################################################
#  File Name:cluster.r
#  Author: xingpengwei
#  Mail: xingwei421@qq.com
#  Created Time: Thu 07 Jan 2021 01:23:41 PM UTC
#################################################

library(pheatmap)
library(edgeR)
library('RColorBrewer')
args <- commandArgs(TRUE)
all_data = read.table(args[1],header=T)
gene_length=all_data[,3]-all_data[,2]
all_data=all_data[,4:length(all_data[1,])]
all_data2=apply(all_data,2,as.numeric)
#data_norm=rpkm(all_data2,gene.length = gene_length)
data_norm=cpm(all_data2,log=FALSE)
data_norm2 = log(data_norm+1,base=2)
a=cor(data_norm2)
pdf("cluster.pdf", width=4.5, height=4)
ann_colors = list(
    FRiP = colorRampPalette(brewer.pal(n =5, name = "BuGn"))(15),
    tss_score = colorRampPalette(brewer.pal(n =5, name = "YlGnBu"))(20)
    )
#pheatmap(a,display_numbers = T,number_color = "black",fontsize_number =9,fontsize_row=15,fontsize_col=15,color = colorRampPalette(brewer.pal(n = 7, name = "RdYlBu"))(100),annotation_row = info,annotation_colors = ann_colors)
pheatmap(a,number_color = "black",color = colorRampPalette(brewer.pal(n = 7, name = "RdBu"))(100),border_color=NA)
dev.off()

#library("FactoMineR")
#library("factoextra")
#res.pca = PCA(t(data_norm))
#pdf(file="PCA.dot.pdf",width=8,height=8)
#fviz_pca_ind(res.pca, pointsize = 5,
# pointshape = 21, fill = "#E7B800",
# repel = TRUE )# Avoid text overlapping (slow if many points)
#dev.off()
save.image('my.Rdata')
