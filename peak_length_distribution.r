#################################################
#  File Name:peak_length_distance.r
#  Author: xingpengwei
#  Mail: xingwei421@qq.com
#  Created Time: Sun 17 Jan 2021 05:32:16 PM UTC
#################################################
library(ggplot2)
library(reshape2)
FACT = read.table("FACT.merge.bed")
normal = read.table("Normal.merge.bed")
encode = read.table("Encode.merge.bed")
FACT_len = FACT$V3-FACT$V2
FACT_len = as.data.frame(FACT_len)
normal_len = normal$V3-normal$V2
normal_len = as.data.frame(normal_len)
encode_len = encode$V3-encode$V2
encode_len = as.data.frame(encode_len)
FACT_MassStat <- hist(FACT_len$FACT_len, breaks=seq(0, ceiling(max(FACT_len)/200)*200, by = 200), plot=FALSE)
FACT_stat = as.data.frame(cbind(FACT_MassStat$breaks,FACT_MassStat$counts))
normal_MassStat <- hist(normal_len$normal_len, breaks=seq(0, ceiling(max(normal_len)/200)*200, by = 200), plot=FALSE)
normal_stat = as.data.frame(cbind(normal_MassStat$breaks,normal_MassStat$counts))
encode_MassStat <- hist(encode_len$encode_len, breaks=seq(0, ceiling(max(encode_len)/200)*200, by = 200), plot=FALSE)
encode_stat = as.data.frame(cbind(encode_MassStat$breaks,encode_MassStat$counts))
normal_stat$V2=normal_stat$V2/sum(normal_stat$V2)
FACT_stat$V2=FACT_stat$V2/sum(FACT_stat$V2)
encode_stat$V2=encode_stat$V2/sum(encode_stat$V2)
all_stat = cbind(FACT_stat[1:500,],normal_stat$V2[1:500],encode_stat$V2[1:500])
all_stat2 = all_stat[1:75,]
colnames(all_stat)=c("length","FACT","normal","encode")
colnames(all_stat2)=c("length","FACT","normal","encode")
all_stat_long = melt(all_stat,id.vars=c("length"))
all_stat_long2 = melt(all_stat2,id.vars=c("length"))
p=ggplot(data = all_stat_long2)+geom_line(aes(x=length,y=value,color=variable))+geom_point(aes(x=length,y=value,color=variable))+theme_bw()+ylab("Percentage of peaks")+xlab("peak length")+scale_x_continuous(breaks = seq(0,15000,1000),expand = c(0,0))+theme(axis.text.x = element_text(angle=45,hjust=0.5,vjust=0.7,color="black"),axis.text.y = element_text(color="black"))
ggsave('all_peak_length_dist.pdf',height=4,width=6)
