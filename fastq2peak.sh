#################################################
#  File Name:fastq2peak.sh
#  Author: xingpengwei
#  Mail: xingwei421@qq.com
#  Created Time: Fri 11 Jun 2021 05:32:16 PM UTC
#################################################
mkdir 01.mapping
mkdir 02.peaks
mkdir 03.bigwig
mkdir 04.TSS_enrich
hs_refanno=/disk1/pengweixing/database/hg38/hg38.refGene.gtf
mm_refanno=/disk1/pengweixing/database/mm10/mm10.refGene.gtf
mm_ref=/disk1/pengweixing/database/mm10/index/mm10.fa
hs_ref=/disk1/pengweixing/database/hg38/index/hg38.fa
currdir=`pwd`
for i in `ls  GM_H3K27ac_ENCODE*gz`
do
out_name=`echo $i|sed s/.fastq.gz//g`
#bowtie2 --end-to-end --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700 -p 20 -x $hs_ref -U $i -S $currdir/01.mapping/${out_name}_bowtie2.sam &> $currdir/01.mapping/${out_name}_bowtie2.txt
#samtools view --threads 10 -S -b $currdir/01.mapping/${out_name}_bowtie2.sam |samtools sort -m 10G - -o $currdir/01.mapping/${out_name}_bowtie2.sort.bam
#/disk1/pengweixing/software/preseq/preseq-3.1.2/build/bin/preseq c_curve -o $currdir/01.mapping/${out_name}.curve.txt -B $currdir/01.mapping/${out_name}_bowtie2.sort.bam
#/disk1/pengweixing/software/preseq/preseq-3.1.2/build/bin/preseq lc_extrap -o $currdir/01.mapping/${out_name}.expect.txt -B $currdir/01.mapping/${out_name}_bowtie2.sort.bam
#java -Xmx4G -jar /home/xingqichen/SOFTWARE/picard-tools-1.119/MarkDuplicates.jar INPUT=$currdir/01.mapping/${out_name}_bowtie2.sort.bam OUTPUT=$currdir/01.mapping/${out_name}_bowtie2.sort.rmdup.bam METRICS_FILE=$currdir/01.mapping/${out_name}_Picard_Metrics_unfiltered_bam.txt VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true REMOVE_DUPLICATES=true &> $currdir/01.mapping/${out_name}.Picard.log
#/disk1/pengweixing/Old_research/snp-calling-pipeline/pipe/sambamba index -t 10 $currdir/01.mapping/${out_name}_bowtie2.sort.bam
#samtools view -h -q 2 -F 0x4 $currdir/01.mapping/${out_name}_bowtie2.sort.bam |awk '$3!="chrM"' |samtools view -b - > $currdir/01.mapping/${out_name}.q2.sort.bam
#samtools flagstat $currdir/01.mapping/${out_name}.q2.sort.bam > $currdir/01.mapping/${out_name}.q2.sort.bam.stat
#sambamba index $currdir/01.mapping/${out_name}.q2.sort.bam
#bamCoverage --numberOfProcessors 20 -b $currdir/01.mapping/${out_name}.q2.sort.bam -o $currdir/03.bigwig/${out_name}.q2.bw   --normalizeUsing CPM
#computeMatrix scale-regions -S  ./03.bigwig/${out_name}.q2.bw \
#                              -R $hs_refanno \
#                              --beforeRegionStartLength 3000 \
#                              --regionBodyLength 2000 \
#                              --afterRegionStartLength 3000 \
#                              --binSize 10 \
#                              --missingDataAsZero \
#                              --sortRegions descend \
#                              --skipZeros -o ./03.bigwig/${out_name}.matrix_gene.mat.gz -p 20
#plotHeatmap --matrixFile ./03.bigwig/${out_name}.matrix_gene.mat.gz --sortRegions no --outFileName ./04.TSS_enrich/${out_name}.scale-regions.enrich.pdf  --colorMap Blues  --heatmapHeight 12 --legendLocation upper-left 
#bedtools bamtobed -i  $currdir/01.mapping/${out_name}.q2.sort.bam >  $currdir/01.mapping/${out_name}.q2.sort.bed
sicer -t $currdir/01.mapping/${out_name}.q2.sort.bed -s hg38 --window_size 200 --gap_size 400 --output_directory 02.peaks
done
