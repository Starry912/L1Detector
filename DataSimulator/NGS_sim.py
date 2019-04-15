#!/usr/bin/env python
# -*- coding: utf-8 -*-

from L1_insertion import *
import pysam
import random
import numpy as np
import sys


#设置参数

ref_path='./resourse/hg19.genome.fa.rm327L1HS.masked.fa'
references=('chr1','chr2')
#references=('chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6')

read_length = 100
isize = 300
read1_file = 'chr_read1.fasta'
read2_file = 'chr_read2.fasta'
l1_prim = 'AGATATACCTAATGCTAGATGACACA'
p_cut=0.8
amplify_max_length=3000
long_amplify_max_times=1
short_amplify_max_times=2
total_base=5*1e8
basic_reads=6.3*1e5*coverage
l1_reads=coverage*total_base/(2*read_length)-basic_reads
if l1_insert_num!=2000:
    l1_amplify_max_times=400
else:
    l1_amplify_max_times=400*1000*l1_reads/(370000.0*l1_insert_num*2)
print 'l1_amplify_max_times',l1_amplify_max_times

def comple(base): #取某个碱基的互补碱基
    if base == 'A' or base == 'a':
        return 'T'
    if base == 'T' or base == 't':
        return 'A'
    if base == 'C' or base == 'c':
        return 'G'
    if base == 'G' or base == 'g':
        return 'C'

def revunit(unit): #将unit逆序
    return unit[::-1]

def compleunit(unit): #取unit的互补链
    # reu = unit[::-1]
    uc = ''
    for i in range(len(unit)):
        uc += str(comple(unit[i]))
    return uc


if __name__ == "__main__":
    #获取参数
    args=sys.argv[1:]
    coverage=float(args[0])
    l1_insert_num=int(args[1])
    #anotation_percent=float(args[2])
    print 'insert_num：',l1_insert_num,' coverage:',coverage
    
    cs_dict=get_cs_dict(cutting_site_file)

    l1_insert_file='anotation_l1_insert' + str(l1_insert_num) + '.txt'
    #anotation_file='anotation_insert.txt'
    l1_insert_dict=get_l1_insert_dict(l1_insert_num,cs_dict,l1_insert_file)

    anotation_file3 = 'anotation_insert3.txt'
    anotation_file5 = 'anotation_insert5.txt'
    anotation_file8 = 'anotation_insert8.txt'

    get_lable_percent_file(l1_insert_file, anotation_file3, 0.3)
    get_lable_percent_file(l1_insert_file, anotation_file5, 0.5)
    get_lable_percent_file(l1_insert_file, anotation_file8, 0.8)
    # 生成read 参数：read_length=100     isize=300

    ref = pysam.FastaFile(ref_path)
    out1=open(read1_file,'w')
    out2=open(read2_file,'w')
    reads_total_num=0
    insert_num=0
    #模拟限制酶剪切 然后得到PE reads
    for chr in references:
        site=sorted(cs_dict[chr])
        head=site[0]
        for i in site:
            if (chr,i) in l1_insert_dict:       #插入L1
                insert_num+=1
                l1_info=l1_insert_dict[(chr,i)]
                l1_seq=l1_info.get_l1_seq()

                cutted_genome_length = int(np.random.normal(2000, 50))
                if cutted_genome_length < isize : cutted_genome_length = isize
                seq_cutted = ref.fetch(chr, i- cutted_genome_length, i).upper()
                seq_inserted= seq_cutted[0:cutted_genome_length- isize] + l1_seq + seq_cutted[cutted_genome_length - isize:] + l1_prim # 切下的序列(l1在中间)
                amplify_times = int(np.random.normal(l1_amplify_max_times, 10))
                #print chr, i,amplify_times
                for n in range(amplify_times):
                        j=random.randint(1,len(seq_inserted))
                        r = int(np.random.normal(0, 5))
                        try:
                            read1 = seq_inserted[r + j: r + j + read_length]
                            read2 = seq_inserted[r + j - read_length + isize: r + j + isize]
                        except:
                            read1 = seq_inserted[len(seq_cutted) - isize: len(seq_cutted) - isize + read_length]
                            read2 = seq_inserted[len(seq_cutted) - read_length:]
                        if len(read1) < read_length * 0.4 or len(read2) < read_length * 0.4: continue
                        out1.write('>' + chr + ':' + str(i) + ':' + str(n) + ':' + str(j) + '\n')  # 写reads
                        out2.write('>' + chr + ':' + str(i) + ':' + str(n) + ':' + str(j) + '\n')
                        out1.write(read1 + '\n')
                        out2.write(compleunit(revunit(read2)) + '\n')
                        reads_total_num += 1
                continue
            iscutted=1 if random.random()< p_cut else 0
            if iscutted==1:
                cutted_genome_length=i-head
                head = i
                if cutted_genome_length < read_length:
                    continue
                if cutted_genome_length > amplify_max_length:
                    continue
                if cutted_genome_length < isize:
                    seq_cutted=ref.fetch(chr, i - cutted_genome_length, i).upper() #+ l1_prim #切下的序列
                    amplify_times=random.randint(0,short_amplify_max_times)
                    for n in range(amplify_times):    #扩增
                        r=random.randint(-5,5)
                        read1 = seq_cutted[r: r + read_length]
                        read2 = seq_cutted[r - read_length:]
                        if len(read1)< read_length * 0.4 or len(read2)< read_length * 0.4 : continue
                        out1.write('>' + chr+':'+ str(i) +':'+ str(n) + '\n')       #写reads
                        out2.write('>' + chr+':'+ str(i) +':'+ str(n) + '\n')
                        out1.write(read1 + '\n')
                        out2.write(compleunit(revunit(read2)) + '\n')
                        reads_total_num+=1
                        continue
                if cutted_genome_length < amplify_max_length:
                    seq_cutted = ref.fetch(chr, i - cutted_genome_length,i).upper() #+ l1_prim  # 切下的序列
                    amplify_times = random.randint(0, long_amplify_max_times)
                    for n in range(amplify_times):
                        for j in range(0,len(seq_cutted),isize):
                            r = random.randint(-5, 5)
                            try:
                                read1 = seq_cutted[r+j: r +j+ read_length]
                                read2 = seq_cutted[r+j - read_length + isize: r +j+ isize]
                            except:
                                read1 = seq_cutted[len(seq_cutted) - isize: len(seq_cutted) - isize + read_length]
                                read2 = seq_cutted[len(seq_cutted) - read_length:]
                            if len(read1) < read_length * 0.4 or len(read2) < read_length * 0.4: continue
                            out1.write('>' + chr + ':' + str(i) + ':' + str(n) + ':' + str(j) + '\n')  # 写reads
                            out2.write('>' + chr + ':' + str(i) + ':' + str(n) + ':' + str(j) + '\n')
                            out1.write(read1+'\n')
                            out2.write(compleunit(revunit(read2))+'\n')
                            reads_total_num += 1
        print chr,':','reads_total_num',reads_total_num,'l1_insert_num',l1_insert_num
    out1.close()
    out2.close()