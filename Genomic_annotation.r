#################################################
#  File Name:SDS_anno2.r
#  Author: xingpengwei
#  Mail: xingwei421@qq.com
#  Created Time: Mon 07 Dec 2020 06:03:34 PM UTC
#################################################
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ChIPseeker)

FACT1 = read.table("GM-T7-ac-1_R1.q2.bed")
colnames(FACT1) = c("chr","start","end")
FACT2 = read.table("GM-T7-ac-2_R1.q2.bed")
colnames(FACT2) = c("chr","start","end")

#THS1 = read.table("THS_Kid_50nM_1.bed")
#colnames(THS1) = c("chr","start","end")
#THS2 = read.table("THS_Kid_50nM_2.bed")
#colnames(THS2) = c("chr","start","end")

ENCODE1 = read.table("GM_H3K27ac_ENCODE_rep1.q2.bed")
colnames(ENCODE1) = c("chr","start","end")
ENCODE2 = read.table("GM_H3K27ac_ENCODE_rep2.q2.bed")
colnames(ENCODE2) = c("chr","start","end")

Normal1 = read.table("GM-normal-ac-1_R1.q2.bed")
colnames(Normal1) = c("chr","start","end")
Normal2 = read.table("GM-normal-ac-2_R1.q2.bed")
colnames(Normal2) = c("chr","start","end")

bbb <- function(data)
{
SDS_data = data
SDS_data=makeGRangesFromDataFrame(SDS_data)
#SDS_data_anno = annotatePeak(SDS_data,tssRegion = c(-3000,3000),TxDb=TxDb.Mmusculus.UCSC.mm9.knownGene)
SDS_data_anno = annotatePeak(SDS_data,tssRegion = c(-3000,3000),TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene)
SDS_anno = SDS_data_anno@annoStat
SDS_anno2 = data.frame(Feature=c("Promoter","5'UTR","3'UTR","Exon","Intron","Dwonstream","Distal Intergenic"),Frequency=c(sum(SDS_anno[1:3,2]),SDS_anno[4,2],SDS_anno[5,2],sum(SDS_anno[6:7,2]),sum(SDS_anno[8:9,2]),SDS_anno[10,2],SDS_anno[11,2]))
SDS_anno2$Feature=factor(SDS_anno2$Feature,levels = c("Promoter","Intron","Distal Intergenic","Dwonstream", "Exon","3'UTR" ,"5'UTR"))
SDS_data_anno@annoStat=SDS_anno2
return(SDS_data_anno)
}

FACT1_anno = bbb(FACT1)
FACT2_anno = bbb(FACT2)

#THS1_anno = bbb(THS1)
#THS2_anno = bbb(THS2)

ENCODE1_anno = bbb(ENCODE1)
ENCODE2_anno = bbb(ENCODE2)

Normal1_anno = bbb(Normal1)
Normal2_anno = bbb(Normal2)

#all = list(FACT1 =FACT1_anno,FACT2 =FACT2_anno, THS1 = THS1_anno, THS2 = THS2_anno, ENCODE1 = ENCODE1_anno, ENCODE2 = ENCODE2_anno,Normal1 = Normal1_anno, Normal2 = Normal2_anno)
all = list(FACT1 =FACT1_anno,FACT2 =FACT2_anno, ENCODE1 = ENCODE1_anno, ENCODE2 = ENCODE2_anno,Normal1 = Normal1_anno, Normal2 = Normal2_anno)

pdf('genomic_Peakanno_all.pdf',width=6,height=3)
plotAnnoBar(all)
dev.off()

