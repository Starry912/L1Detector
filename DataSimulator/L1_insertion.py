#!/usr/bin/env python
# -*- coding: utf-8 -*-

from L1_creater import *
import pysam
import random


#用于生成fasta

#参考的剪切位点文件
cutting_site_file='./resourse/vectorette-enzyme-cutting-site-info.txt.wgsrh'
ref_path='./resourse/hg19.genome.fa.rm327L1HS.masked.fa'

#指定插入chr
references=('chr1', 'chr2')
# references=('chrM', 'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY')
#l1_insert_num=2000


class cutting_site_info:
    def set_cs(self,list):
        self.chr = list[0].split('-')[0]
        self.seq = list[1].split('-')[1]
        self.type = list[2]
        self.hits = list[3].split('-')[1]
        return self.chr



def get_cs_dict(cutting_site_file):
    # 读取cs
    cs_dict = {}
    cs = cutting_site_info()
    with open(cutting_site_file, 'r') as f1:
        for line in f1:
            if line[0] == '>':
                chr = cs.set_cs(line[1:].split(':'))
                if chr in cs_dict.keys():
                    continue
                else:
                    cs_dict[chr] = []
                    continue
            try:
                cs_dict[chr] = cs_dict[chr] + map(int, line.strip('\n').split(','))
            except:
                pass
    return cs_dict


class l1_insert_info:
    def set_l1(self,list):
        self.chr = list[0]
        self.target_site = list[1]
        self.type = list[2]
        self.l1_seq = list[3]

        #self.TSD = list[5]

    def get_l1_seq(self):
        return self.l1_seq



def get_l1_insert_dict(l1_insert_num,cs_dict,l1_insert_file):
    #l1插入位点list
    l1_insert_dict={}
    tanget_site_list={}
    l1_info = l1_insert_info()
    fout1=open(l1_insert_file, 'w')
    #fout2=open(anotation_file_path,'w')
    each_l1_insert_num=l1_insert_num/2
    for chr in references:
        tanget_site_list[chr] = random.sample(cs_dict[chr],each_l1_insert_num)
    for i in range(each_l1_insert_num):
        for chr in references:
            k=random.choice([1,2,2,2,3])      #类型
            #chr=random.choice(references)
            tanget_site=tanget_site_list[chr][i]
            tanget_site_start=tanget_site-300
            site=(chr,tanget_site)
            l1_seq = l1_creater(k, chr, tanget_site_start, ref_path, l1_file_path)
            l1_info.set_l1([chr,tanget_site,k,l1_seq])

            l1_insert_dict[site]=l1_info
            l1_info = l1_insert_info()
            print >>fout1,chr,tanget_site_start,tanget_site,k
            #if i%10<anotation_percent*10:
            #    print >> fout2, chr, tanget_site_start, tanget_site, k
    fout1.close()
    #fout2.close()
    return l1_insert_dict


def get_lable_percent_file(l1_insert_file,anotation_file,anotation_percent):
    """获得label"""
    i=0
    with open(l1_insert_file,'r') as f:
        with open(anotation_file,'w') as f_out:
            for line in f:
                if i%10<anotation_percent*10:
                    f_out.write(line)
                i+=1