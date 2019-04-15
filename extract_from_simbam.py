#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import pysam
import numpy as np
import re
from math import *
from datetime import datetime
from read_function import *
#from PCA import *
#from LVQ import *
from fliter import *
#from semi_Kmeans import *


l1_insert_num=int(sys.argv[1])
window_size=int(sys.argv[2])
anotation_file_path=sys.argv[3]
print 'window_size',window_size

human_bam_path='./data/0102_hg_sort.bam'
#human_bam_path = '/mnt/hgfs/transopon_data/data/output/bowtie_human/A152-Normal.s5_9_Sl31275..fastq.cleaned.fastq.pcsort.bam'
l1hs_bam_path = './data/0102_l1.sorted.bam'
log_path = '/mnt/hgfs/transposon/'+datetime.now().strftime('%Y-%m-%d %H:%M:%S')+'.txt'
#anotation_file_path = '/mnt/hgfs/transopon_data/annotation/ucsc.rm327.rmsk.l1HS.coord.with-strand.txt.Ta.only.bed'
#anotation_file_path = 'anotation_insert.txt'
postive_gnd_file='anotation_l1_insert'+str(l1_insert_num)+'.txt'


prim='AGATATACCTAATGCTAGATGACACA'
i_size=300
i_size_std=30


f_out=open(log_path,'w')
# step1 提取candidate bk
print "开始step1"
print  datetime.now().strftime('%Y-%m-%d %H:%M:%S')

f_human = pysam.AlignmentFile(human_bam_path, 'rb')

# bamfile.nreferences /bamfile.lengths /bamfile.references
candidate_area_dict = {}
chr_count = 0
window_count=0

for chr in f_human.references[0:6]:
    if chr == 'chrM':
        chr_count += 1
        continue
    for length in range(0, f_human.lengths[chr_count], window_size):
        candidate_read_count = 0
        if length + window_size < f_human.lengths[chr_count]:
            iter_human =  f_human.fetch(chr, length, length + window_size)
            for r in iter_human:
                if not r.is_unmapped:
                    # 是否是J read
                    if re.match('[0-9]+S(.*)[0-9]+M', r.cigarstring) or re.match('[0-9]+M(.*)[0-9]+S', r.cigarstring):
                        # r_l1hs_pairs = f_l1hs_index.find(r.qname)
                        # for r_l1hs in r_l1hs_pairs:
                        #     if r_l1hs.flag == r.flag:
                        #         if re.match('[0-9]+S(.*)[0-9]+M', r_l1hs.cigarstring) or re.match('[0-9]+M(.*)[0-9]+S',
                        #                                                                       r_l1hs.cigarstring):
                        candidate_read_count+=1
            if candidate_read_count!=0:
                # 计算（start,end)
                (start, end) = getMaxBlocks(f_human.fetch(chr, length, length + window_size))
                candidate_area_dict[chr, start, end] = candidate_read_count
                window_count += 1
        else:
            iter_human = f_human.fetch(chr, length, f_human.lengths[chr_count])
            for r in iter_human:
                if not r.is_unmapped:
                    # 是否是J read
                    if re.match('[0-9]+S(.*)[0-9]+M', r.cigarstring) or re.match('[0-9]+M(.*)[0-9]+S', r.cigarstring):
                        # r_l1hs_pairs = f_l1hs_index.find(r.qname)
                        # for r_l1hs in r_l1hs_pairs:
                        #     if r_l1hs.flag == r.flag:
                        #         if re.match('[0-9]+S(.*)[0-9]+M', r_l1hs.cigarstring) or re.match('[0-9]+M(.*)[0-9]+S',
                        #                                                                       r_l1hs.cigarstring):
                        candidate_read_count += 1
            if candidate_read_count != 0:
                #计算（start,end)
                (start, end) = getMaxBlocks(f_human.fetch(chr, length, length + window_size))
                candidate_area_dict[chr, start, end] = candidate_read_count
                window_count += 1
        # if len(candidate_area_dict) > 2000:
        #     print candidate_area_dict
        #     break
    print chr + "finish",window_count
    chr_count += 1
    # if chr_count > 1:
    #     break


print "完成step1"
print  datetime.now().strftime('%Y-%m-%d %H:%M:%S')

# step2 filter
adjust(candidate_area_dict,f_human)
merge(candidate_area_dict,window_size)

# step3 提取特征 标签
# width /depth /variant index /pA /J
# flag


N = len(candidate_area_dict)
print "candidate site",N

# f_l1hs = pysam.AlignmentFile(l1hs_bam_path, 'rb')
# f_l1hs_index = pysam.IndexedReads(f_l1hs)
# f_l1hs_index.build()
# print >> f_out, "f_l1hs_index.build"
pos=[]
feature = np.zeros([N, 13])
label = np.zeros([N, 1])
gnd = np.zeros([N, 1])
postive_list=AnnotationPostive(anotation_file_path)
postive_gnd_list=AnnotationPostive(postive_gnd_file)
count = 0
for i in candidate_area_dict:
    chr = i[0]
    start = i[1]
    end = i[2]
    pos.append(i)
    #feature
    prim_count = 0
    #varition
    #feature3 = sum(r.get_tag('NM') for r in f_human.fetch(chr, start, stop) if r.has_tag('NM'))/f_human.count(chr, start, stop)
    #polyA
    #feature4 = ployARate(candidate_area_dict[i]) #+ployTRate(candidate_area_dict[i])
    #N split
    #feature5 = len(candidate_area_dict[i])
    #
    width = end - start
    total_read = f_human.count(chr, start, end)
    coverage = np.sum(f_human.count_coverage(chr, start, end)) / float(end - start)
    ployA_read=0
    low_mapping_quality_read = 0
    varition_read=0  # %m%s
    abnormal_insert_size_read=0
    split_mapped_read = 0
    incomplete_mapped_read = 0
    part_map_to_l1_read = 0
    mate_map_to_other_read = 0
    mate_map_to_l1_read = 0
    abnormal_stand_read = 0


    for r in f_human.fetch(chr, start, end):
        if not r.is_unmapped:
            if not r.mate_is_unmapped:
                if r.isize > 320:
                    abnormal_insert_size_read += 1
                if r.isize+f_human.mate(r).isize!=0:
                    abnormal_stand_read+=1
                if r.mapping_quality < 20:
                    low_mapping_quality_read +=1
                if r.reference_name != f_human.mate(r).reference_name:
                    mate_map_to_other_read += 1
                    # for r_l in f_l1hs_index.find(f_human.mate(r).qname):
                    #     if not r_l.is_unmapped:
                    #         mate_map_to_l1_read += 1
                if re.match('[0-9]+S(.*)[0-9]+M', r.cigarstring) or re.match('[0-9]+M(.*)[0-9]+S', r.cigarstring):
                    incomplete_mapped_read+=1
                    if ('AATAAA' in r.seq) or ('ATTAAA' in r.seq) or ('TTATTT' in r.seq) or ('TAATTT' in r.seq):
                        ployA_read += 1
                    consensus, consensus_len = findPrim(prim, r.seq)
                    if consensus_len > 0.6 * len(prim):
                        prim_count += 1
                if re.match('(.*)[0-9]+M(.*)[0-9]+M(.*)', r.cigarstring):
                    split_mapped_read+=1
                # for r_l in f_l1hs_index.find(r.qname):
                #     if not r_l.is_unmapped:
                #         part_map_to_l1_read+=1

        if r.has_tag('NM') and r.get_tag('NM')!=0:
            varition_read+=1



    if ployA_read==0:
        ployA_rate=0
    else:
        ployA_rate=ployA_read/float(incomplete_mapped_read)
    #feature[count] =[log(feature1,2) if feature1>0 else feature1, log(feature2,2) if feature2>0 else feature2, log(feature3,2) if feature3>0 else feature3, feature4, log(feature5,2)]
    feature[count] = [width,total_read,coverage,ployA_rate,low_mapping_quality_read, varition_read, abnormal_insert_size_read,
                      split_mapped_read,incomplete_mapped_read,part_map_to_l1_read, mate_map_to_other_read, mate_map_to_l1_read,
                      abnormal_stand_read ]


    [label[count], p] = FindType(count, i, postive_list, prim_count)
    if p != 0: label[p] = -1
    [gnd[count], p] = CheckPos(count, i, postive_gnd_list)
    if p != 0: label[p] = -1
    #计算flag
    # prim_count=0
    # for r in candidate_area_dict[i]:
    #     consensus,consensus_len =findPrim(prim,r.seq)
    #     if consensus_len > 0.6*len(prim):
    #         prim_count+=1
    # prim_count=sum(1 for r in candidate_area_dict[i] if findPrim(prim,r.seq)[1] > 10)
    # [label[count], p] = FindType(count, i, postive_list, prim_count)
    #
    count += 1
    if count%100==0: print count

print "完成step2"
print  datetime.now().strftime('%Y-%m-%d %H:%M:%S')

f_human.close()
#f_l1hs.close()
np.savez('feature.npz',pos,feature,label,gnd)
# write_result(pos,feature,label,gnd)
#
#
# U=feature[np.where(label==0)[0]]
# L=np.hstack((feature[np.where(label!=0)[0]],label[np.where(label!=0)[0]]))


#cluster=semi_kMeans(L, U, distMeas=distEclud, initial_centriod=newCent)
# np.save('cluster.npy',cluster)

#print accuracy(gnd,cluster)

# excel = xlwt.Workbook(encoding='utf-8', style_compression=0)
# sheet = excel.add_sheet('result', cell_overwrite_ok=True)
# for i in range(len(cluster)):
#     sheet.write(i, 0, pos[i][0])
#     sheet.write(i, 1, pos[i][1])
#     sheet.write(i, 2, pos[i][2])
#     sheet.write(i, 3, int(gnd[i]))
#     sheet.write(i, 4, int(cluster[i]))
#
# excel.save(r'/mnt/hgfs/transposon/result.xls')
# step4 PCA
    #1 标准化 ，定标 (Scale, 数值除以标准差) 中心化（centering，数值减去平均）
# feature_scaled=np.divide(feature,np.std(feature,axis=0))
# feature_centered=np.subtract(feature_scaled,np.mean(feature,axis=0))
    #2 PCA
#lowDDataMat, reconMat = pca(feature_centered,0.9)
# step5 lVQ
# P = LVQ(L,feature,label, 0.01, 5000)
# l2 = pridect(np.asarray(f), T, np.asarray(P))
# accuracy(l,l2)
# step6 pridect
#pridect(np.asarray(lowDDataMat),np.asarray(P))

f_out.close()


