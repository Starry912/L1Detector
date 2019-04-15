#!/usr/bin/env python
# -*- coding: utf-8 -*-

from read_function import *
import numpy as np

from sklearn.metrics import classification_report

l1_insert_num=2000
postive_gnd_file='./anotation_l1_insert'+str(l1_insert_num)+'.txt'
postive_gnd_list=AnnotationPostive3(postive_gnd_file)

result_file='/mnt/hgfs/transopon_data/data/out/model/chr1_read..fasta.cleaned.fasta.pcsort.qyname.bam.L1HSAligned.sc.leftvsright.bed.consensusbp.wsize100.regwsize1.minreads1.clip1.clipflk5.mindis150.rmskta.uniqgs.bed.csinfo.lm.l1hs.forsf'
#
result_list=[]
with open(result_file,'r') as result:
    for line in result:
        if line[0]=='H': continue
        TIP=line.strip('\n').split('\t')
        result_list.append([TIP[0], int(TIP[5]), int(TIP[6]),1])

p=np.zeros(len(result_list))
g=np.zeros(len(result_list))
count=0
for a in result_list:
    p[count]=a[3]
    [g[count],m]=CheckPos(count, a, postive_gnd_list)
    if m != 0: g[m] = 0
    count+=1


for b in postive_gnd_list:
    if b[3]==0:
        p=np.hstack((p,0))
        g = np.hstack((g, b[5]))

# target_names = ['nagetive','positive']
# print(classification_report(g,p,target_names=target_names))
print 'final accuracy: %.2f' % (np.mean(g == p))
accuracy_score(g,p)
recall_score(g,p)
confusion_matrix(g,p).ravel()



result_file='./data/group1_l1_insert_num/result(2).txt'
result_list=[]
with open(result_file,'r') as result:
    for line in result:
        TIP=line.strip('\n').split(' ')
        result_list.append([TIP[0], int(TIP[1]), int(TIP[2]), int(TIP[-1])])
        # 0 chr 1 start 2 end 3 type

p=np.zeros(len(result_list))
g=np.zeros(len(result_list))
count=0
for a in result_list:
    p[count]=a[3]
    [g[count],m]=CheckPos(count, a, postive_gnd_list)
    if m != 0: g[m] = 0
    count+=1


for b in postive_gnd_list:
    if b[3]==0:
        p=np.hstack((p,0))
        g = np.hstack((g, b[5]))
print 'final accuracy: %.2f' % (np.mean(g == p))
target_names = ['nagetive', 'full-length', 'trunctive','conductive']
print(classification_report(g,p,target_names=target_names))
print 'final accuracy: %.2f' % (np.mean(g == p))


