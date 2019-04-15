#!/usr/bin/env python
# -*- coding: utf-8 -*-

import random
import pysam

#ref and L1 file path
ref_path='./resourse/hg19.genome.fa.rm327L1HS.masked.fa'
l1_file_path='./resourse/Homo_sapiens_L1.L1HS.fa'
l1_length=6064



def l1_creater(k,chr,tanget_site,ref_path,l1_file_path):
    """生成一个l1插入单元   class1: 5'-truncation   class2:5'truncation with an inversion class3:3'transduction"""
    ref=pysam.FastaFile(ref_path)
    l1=pysam.FastaFile(l1_file_path)
    TSD_length=random.choice(range(9,12))
    TSD=ref.fetch(chr,tanget_site,TSD_length+tanget_site).upper()
    ORF2_length=random.choice(xrange(300,2000))
    if k==1:
        return l1_class1_creater(TSD,l1)
    if k==2:
        return l1_class2_creater(TSD,ORF2_length,l1)
    # if k==3:
    #     ORF1_length=random.choice(xrange(100,500))
    #     return l1_class3_creater(TSD,ORF1_length,ORF2_length,l1)
    if k==3:
        P_chr=random.choice(ref.references)
        P_TSD=ref.fetch(P_chr,tanget_site-TSD_length,tanget_site)
        P_DNA_length=random.choice(xrange(100,400))
        P_DNA_site=random.choice(xrange(1000,10000000))
        P_DNA=ref.fetch(P_chr,P_DNA_site,P_DNA_site+P_DNA_length).upper()
        pA='AATAAAATTAAA'
        return l1_class4_creater(TSD,ORF2_length,P_TSD,P_DNA,pA,l1)


def l1_class1_creater(TSD,l1):
    """class1  TSD+l1+TSD"""
    l1_full=l1.fetch('L1HS',0,l1_length)
    l1_class1=TSD+l1_full #+TSD
    return l1_class1.upper()


def l1_class2_creater(TSD,ORF2_length,l1):
    """class1  TSD+ORF2+3'UTR+TSD"""
    ORF2=l1.fetch('L1HS',l1_length-ORF2_length,l1_length)
    l1_class2=TSD+ORF2 #+TSD
    return l1_class2.upper()


def l1_class3_creater(TSD,ORF1_length,ORF2_length,l1):
    """class2 TSD+//RO+//F2+3'UTR+TSD"""
    ORF1_reverse=l1.fetch('L1HS',l1_length - ORF2_length - ORF1_length,l1_length - ORF2_length)
    ORF2 = l1.fetch('L1HS', l1_length - ORF2_length, l1_length)
    l1_class3=TSD+ORF1_reverse+ORF2 #+TSD
    return l1_class3.upper()


def l1_class4_creater(TSD,ORF2_length,P_TSD,P_DNA,pA,l1):
    """class3 TSD+ORF2+3'UTR+P_TSD+3'P_DNA+pA+TSD"""
    ORF2 = l1.fetch('L1HS', l1_length - ORF2_length, l1_length)
    l1_class4=TSD+ORF2+P_TSD+P_DNA+pA #+TSD
    return l1_class4.upper()

if __name__ == "__main__":
    l1_1=l1_creater(1,'chr1',1000,ref_path,l1_file_path)
    l1_2 = l1_creater(2, 'chr1', 1000, ref_path, l1_file_path)
    l1_3 = l1_creater(3, 'chr1', 1000, ref_path, l1_file_path)
    print l1_1
    print l1_2
    print l1_3
