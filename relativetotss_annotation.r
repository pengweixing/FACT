#################################################
#  File Name:relativetotss.r
#  Author: xingpengwei
#  Mail: xingwei421@qq.com
#  Created Time: Sat 23 Jan 2021 02:58:06 PM UTC
#################################################

library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
args = commandArgs(T)
#txdb <- TxDb.Mmusculus.UCSC.mm9.knownGene
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
#files=list(T7="KM.T7.merge.bed")
#files=list(FACT1='FFPE_kidney_1.bed',FACT2='FFPE_kidney_2.bed',THS1='THS_Kid_50nM_1.bed',THS2='THS_Kid_50nM_2.bed',ENCODE1='Kidney_SD_1.bed',ENCODE2='Kidney_SD_2.bed',Nomal1='NAK_1.bed',Nomal2='NAK_2.bed')
files=list(FACT1='GM-T7-ac-1_R1.q2.bed',FACT2='GM-T7-ac-2_R1.q2.bed', ENCODE1='GM_H3K27ac_ENCODE_rep1.q2.bed',ENCODE2='GM_H3K27ac_ENCODE_rep2.q2.bed',Nomal1='GM-normal-ac-1_R1.q2.bed',Nomal2='GM-normal-ac-2_R1.q2.bed')
peakAnnoList <- lapply(files, annotatePeak, TxDb=txdb, tssRegion=c(-3000, 3000), verbose=FALSE)
pdf(file=paste0("relative_to_tss.pdf"),height= 3, width = 6)
plotDistToTSS(peakAnnoList)
dev.off()

