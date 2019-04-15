#!/usr/bin/env bash

l1_insert_num=(1000 2000 5000 10000)
label_percent=(3 5 8)  #0.3 0.5 0.8
coverage=(0.5 0.8 1 1.2)
window_size=(200 500 1000)

#

#group1_l1_insert_num
#mkdir ./data/group1_l1_insert_num
#for ((group=0;group<4;group++));
#do
#echo ${l1_insert_num[$group]}
#
#python ./simulation/NGS_sim.py 0.8 ${l1_insert_num[$group]}
#
#bowtie2 --very-fast-local -f -x   /mnt/hgfs/transopon_data/data/reference_seqshg19_rmsk-L1HS_masked/hg19.genome.fa.rm327L1HS.masked -1 ./chr_read1.fasta -2 ./chr_read2.fasta -S ./data/0102_hg.sam
#
#samtools view -b -S -t /mnt/hgfs/transopon_data/data/reference_seqshg19_rmsk-L1HS_masked/hg19.genome.fa.rm327L1HS.masked.fa.fai ./data/0102_hg.sam > ./data/0102_hg.bam
#samtools sort -l 9 -o ./data/0102_hg_sort.bam ./data/0102_hg.bam
#samtools index ./data/0102_hg_sort.bam
#
#python ./ExtractFromSimBam.py ${l1_insert_num[$group]} 500 anotation_insert5.txt
#python ./Model.py
##bowtie2 -f -x  /mnt/hgfs/transopon_data/data/reference_seqsL1HS_RepBase/Homo_sapiens_L1.L1HS  -1 ./chr_read1.fasta -2 ./chr_read2.fasta -S ./data/0102_l1.sam
#
##samtools view -b -S -t /mnt/hgfs/transopon_data/data/reference_seqsL1HS_RepBase/Homo_sapiens_L1.L1HS.fa.fai ./data/0102_l1.sam > ./data/0102_l1.bam
##samtools sort -l 9 -o ./data/0102_l1.sorted.bam ./data/0102_l1.bam
##samtools index ./data/0102_l1.sorted.bam
#zip ./data/group1_l1_insert_num/group1v$group.zip ./data/0102_hg_sort.bam ./anotation_l1_insert${l1_insert_num[$group]}.txt ./result.txt ./metrics.txt
##zip ./data/group1_l1_insert_num/group1v$group.zip ./data/0102_hg_sort.bam ./anotation_l1_insert1000.txt ./result.txt ./metrics.txt
#done

#group2_label_percent
#mkdir ./data/group2_label_percent
#
#python ./simulation/NGS_sim.py 0.8 2000
#
#bowtie2 --very-fast-local -f -x   /mnt/hgfs/transopon_data/data/reference_seqshg19_rmsk-L1HS_masked/hg19.genome.fa.rm327L1HS.masked -1 ./chr_read1.fasta -2 ./chr_read2.fasta -S ./data/0102_hg.sam
#
#samtools view -b -S -t /mnt/hgfs/transopon_data/data/reference_seqshg19_rmsk-L1HS_masked/hg19.genome.fa.rm327L1HS.masked.fa.fai ./data/0102_hg.sam > ./data/0102_hg.bam
#samtools sort -l 9 -o ./data/0102_hg_sort.bam ./data/0102_hg.bam
#samtools index ./data/0102_hg_sort.bam
#
#for ((group=0;group<3;group++));
#do
#echo ${label_percent[$group]}
#
#python ./ExtractFromSimBam.py 2000 500 anotation_insert${label_percent[$group]}.txt
#python ./Model.py
#
#zip ./data/group2_label_percent/group2v$group.zip ./data/0102_hg_sort.bam ./anotation_l1_insert2000.txt ./result.txt ./metrics.txt ./result.npz ./feature.npz
#done
#
##group3_coverage
#mkdir ./data/group3_coverage
#for ((group=0;group<4;group++));
#do
#echo ${coverage[$group]}
#python ./simulation/NGS_sim.py ${coverage[$group]} 2000
#
#bowtie2 --very-fast-local -f -x   /mnt/hgfs/transopon_data/data/reference_seqshg19_rmsk-L1HS_masked/hg19.genome.fa.rm327L1HS.masked -1 ./chr_read1.fasta -2 ./chr_read2.fasta -S ./data/0102_hg.sam
#
#samtools view -b -S -t /mnt/hgfs/transopon_data/data/reference_seqshg19_rmsk-L1HS_masked/hg19.genome.fa.rm327L1HS.masked.fa.fai ./data/0102_hg.sam > ./data/0102_hg.bam
#samtools sort -l 9 -o ./data/0102_hg_sort.bam ./data/0102_hg.bam
#samtools index ./data/0102_hg_sort.bam
#
#python ./ExtractFromSimBam.py 2000 500 anotation_insert5.txt
#python ./Model.py
#
#zip ./data/group3_coverage/group3v$group.zip ./data/0102_hg_sort.bam ./anotation_l1_insert2000.txt ./result.txt ./metrics.txt ./result.npz ./feature.npz
#done

##group4_window_size
#mkdir ./data/group4_window_size
#python ./simulation/NGS_sim.py 0.5 2000
#
#bowtie2 --very-fast-local -f -x   /mnt/hgfs/transopon_data/data/reference_seqshg19_rmsk-L1HS_masked/hg19.genome.fa.rm327L1HS.masked -1 ./chr_read1.fasta -2 ./chr_read2.fasta -S ./data/0102_hg.sam
#
#samtools view -b -S -t /mnt/hgfs/transopon_data/data/reference_seqshg19_rmsk-L1HS_masked/hg19.genome.fa.rm327L1HS.masked.fa.fai ./data/0102_hg.sam > ./data/0102_hg.bam
#samtools sort -l 9 -o ./data/0102_hg_sort.bam ./data/0102_hg.bam
#samtools index ./data/0102_hg_sort.bam
#for ((group=0;group<3;group++));
#do
#echo ${window_size[$group]}
#
#python ./ExtractFromSimBam.py 2000 ${window_size[$group]} anotation_insert5.txt
#python ./Model.py
#
#zip ./data/group4_window_size/group4v$group.zip ./data/0102_hg_sort.bam ./anotation_l1_insert2000.txt ./result.txt ./metrics.txt ./result.npz ./feature.npz
#done

#gorup5_constrast1
#mkdir ./data/gorup5_constrast1
#python ./simulation/NGS_sim.py 0.8 2000

#my
#bowtie2 --local --sensitive -f -x   /mnt/hgfs/transopon_data/data/reference_seqshg19_rmsk-L1HS_masked/hg19.genome.fa.rm327L1HS.masked -1 ./chr_read1.fasta -2 ./chr_read2.fasta -S ./data/0102_hg.sam
#samtools view -b -S -t /mnt/hgfs/transopon_data/data/reference_seqshg19_rmsk-L1HS_masked/hg19.genome.fa.rm327L1HS.masked.fa.fai ./data/0102_hg.sam > ./data/0102_hg.bam
#samtools sort -l 9 -o ./data/0102_hg_sort.bam ./data/0102_hg.bam
#samtools index ./data/0102_hg_sort.bam
#python ./ExtractFromSimBam.py 2000 300 anotation_insert5.txt
#zip ./gorup5_constrast1/group5v$group.zip ./data/0102_hg_sort.bam ./anotation_l1_insert2000.txt ./result.txt ./metrics.txt

#L1seqhunter
#./TIPseqHunterPipelineJar.sh ./data/fqfiles/ ./data/output/ A152-Normal.s5_9_SL31275.R1.fastq R1 R2 70197523 outnamesortbamfile=A152-Normal.s5_9_SL31275..fastq.cleaned.fastq.pcsort.qyname.bam
/mnt/hgfs/transopon_data/TIPseqHunterPipelineJar.sh /mnt/hgfs/transopon_data/data/fq /mnt/hgfs/transopon_data/data/out  chr1_read.R1.fasta R1 R2 1386412
