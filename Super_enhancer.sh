#################################################
#  File Name:super2.sh
#  Author: xingpengwei
#  Mail: xingwei421@qq.com
#  Created Time: Sat 20 Feb 2021 10:58:39 AM UTC
#################################################

out=HGBM_p2
se=HG38
bedtools bamtobed -i  FA-HG-A5000.merge.nochrM.bam| awk 'BEGIN{{OFS="\t"}}{{$4="N";$5="1000";print $0}}' > $out.bed
macs2 callpeak -t $out.bed -f BED -n $out -g hs --nomodel --shift 0 --extsize 150  --keep-dup all -B --SPMR  -p 0.00001
super enhancer calling###
python /disk1/pengweixing/software/rose/bed2gff.py ${out}_peaks.narrowPeak ${out}.peak.gff
samtools index FA-HG-A5000.merge.nochrM.bam
/usr/bin/python /disk1/pengweixing/software/rose/ROSE_main.py -i ${out}.peak.gff -r FA-HG-A5000.merge.nochrM.bam -o $out.enhancer -g $se &
#/usr/bin/python /disk1/pengweixing/software/rose/ROSE_main.py -i merge.peak.clean.gff -r merge.bam -o merge.peak.clean -g $se
