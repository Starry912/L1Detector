#!/usr/bin/env python
# -*- coding: utf-8 -*-

import re
import xlwt



def getMaxBlocks(read_list):
    start=0
    end=0
    for r in read_list:
        if start == 0:
            start=r.reference_start
            end=r.aend
        else:
            start=min(r.reference_start,start)
            end=max(r.aend,end)
    return (start,end)

def findPrim(prim,seq):
    #序列中是否包含引物序列
    consensus=''
    for i in range(len(prim)):
        for j in range(i+1,len(prim)+1):
            temp=prim[i:j]
            if seq.find(temp)<0:
                break
            elif len(consensus)<len(temp):
                consensus=temp
    return [consensus,len(consensus)]

def hamming(seq1,seq2):
    mutation = [i for i in range(len(seq1)) if seq1[i] != seq2[i]]
    return len(mutation)

def matchCount(cigarstring):
    return sum(map(int,re.findall('([0-9]+)M',cigarstring)))

# 计算ployA%
def ployARate(r_list):
    totol_ployA = 0
    ployA = 0
    for r in r_list:
        if re.match('[0-9]+S(.*)[0-9]+M', r.cigarstring):  # is 3' read  [0-9]+S(.*)[0-9]+M
            totol_ployA += 1
            if ('AATAAA' in r.seq) or ('ATTAAA' in r.seq):
                ployA += 1
    if totol_ployA == 0:
        return 0
    else:
        return ployA/float(totol_ployA)

# 计算ployT%
def ployTRate(r_list):
    totol_ployT = 0
    ployT = 0
    for r in r_list:
        if re.match('[0-9]+S(.*)[0-9]+M', r.cigarstring):  # is 3' read  [0-9]+S(.*)[0-9]+M
            totol_ployT += 1
            if ('TTATTT' in r.seq) or ('TAATTT' in r.seq):
                ployT += 1
    if totol_ployT == 0:
        return 0
    else:
        return ployT/float(totol_ployT)


def mismathCount(r_list):
    # 计算mismatch
    mismatch = 0
    for r in r_list:
        MD = r.get_tag('MD')
        pat = "[0-9]+[ATGC]+"
        MD_list = re.findall(pat, MD)
        for m in MD_list:
            for n in m:
                if n == 'A' or n == 'T' or n == 'G' or n == 'C':
                    mismatch += 1
    return mismatch

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

def AnnotationPostive(anotation_file_path):
    postive_list = []
    with open(anotation_file_path, 'r') as f:
        for line in f:
            TIP = line.strip().split(' ')
            postive_list.append([TIP[0], int(TIP[1]), int(TIP[2]),0,0,int(TIP[3])])
            #0 chr 1 start 2 end 3 pos in candidate  4 distance 5 type
    return postive_list

def AnnotationPostive2(anotation_file_path):
    postive_list = []
    with open(anotation_file_path, 'r') as f:
        for line in f:
            TIP = line.strip().split('\t')
            postive_list.append([TIP[0], int(TIP[1]), int(TIP[2]),0,0,TIP[3]])
            #0 chr 1 start 2 end 3 pos in candidate  4 distance 5 type
    return postive_list

def AnnotationPostive3(anotation_file_path):
    postive_list = []
    with open(anotation_file_path, 'r') as f:
        for line in f:
            TIP = line.strip().split(' ')
            postive_list.append([TIP[0], int(TIP[1]), int(TIP[2]),0,0,1])
            #0 chr 1 start 2 end 3 pos in candidate  4 distance 5 type
    return postive_list


# def FindType(count,a,postive_list,prim_count):
#     """"ascertain lable"""
#     #if prim_count==0: return [0,0]
#     for b in postive_list:
#         if b[0]==a[0]:
#             distance = abs(a[1] - b[1]) + abs(a[2] - b[2])
#             width=a[2]-a[1]
#             if distance < 0.0001 * b[1]:
#                 n=postive_list.index(b)
#                 if postive_list[n][3]==0:
#                     postive_list[n][3] = count
#                     postive_list[n][4] = width
#                     return [1,0]
#                 if postive_list[n][3]!=0 and postive_list[n][4]<width:
#                     p=postive_list[n][3]
#                     postive_list[n][3]=count
#                     postive_list[n][4]=width
#                     return [1,p]
#     if prim_count == 0: return [0, 0]
#     return [2,0]


def FindType(count,a,postive_list,prim_count):
    """"ascertain lable"""
    #if prim_count==0: return [0,0]
    for b in postive_list:
        if b[0]==a[0]:
            distance = (abs(a[1] - b[1]) + abs(a[2] - b[2]))/2
            width=a[2]-a[1]
            if distance < 3 * width:
                n=postive_list.index(b)
                if postive_list[n][4]==0:
                    postive_list[n][3] = count
                    postive_list[n][4] = distance
                    return [postive_list[n][5],0]
                if postive_list[n][4]!=0 and postive_list[n][4]<distance:
                    p=postive_list[n][3]
                    postive_list[n][3]=count
                    postive_list[n][4]=distance
                    return [postive_list[n][5],p]
    if prim_count == 0: return [0, 0]
    return [-1,0]

def FindType2(count,a,postive_list,prim_count):
    """"ascertain lable"""
    #if prim_count==0: return [0,0]
    for b in postive_list:
        if b[0]==a[0]:
            distance = (abs(a[1] - b[1]) + abs(a[2] - b[2]))/2
            width=a[2]-a[1]
            if distance < 3 * width:
                n=postive_list.index(b)
                if postive_list[n][4]==0:
                    postive_list[n][3] = count
                    postive_list[n][4] = distance
                    return [1,0]
                if postive_list[n][4]!=0 and postive_list[n][4]<distance:
                    p=postive_list[n][3]
                    postive_list[n][3]=count
                    postive_list[n][4]=distance
                    return [1,p]
    if prim_count == 0: return [0, 0]
    return [-1,0]


def CheckPos(count,a,postive_list):
    """"ascertain lable"""
    #if prim_count==0: return [0,0]
    for b in postive_list:
        if b[0]==a[0]:
            distance = (abs(a[1] - b[1]) + abs(a[2] - b[2]))/2
            width=a[2]-a[1]
            if distance < 3 * width:
                n=postive_list.index(b)
                if postive_list[n][4]==0:
                    postive_list[n][3] = count
                    postive_list[n][4] = distance
                    return [postive_list[n][5],0]
                if postive_list[n][4]!=0 and postive_list[n][4]<distance:
                    p=postive_list[n][3]
                    postive_list[n][3]=count
                    postive_list[n][4]=distance
                    return [postive_list[n][5],p]
    return [0,0]




def write_result(pos,feature,label,gnd):
    excel = xlwt.Workbook(encoding='utf-8', style_compression=0)
    sheet=excel.add_sheet('myfeature',cell_overwrite_ok=True)
    for i in range(len(feature)):
        for j in range(len(pos[i])):
            sheet.write(i,j, pos[i][j])
        for j in range(len(feature[i])):
            sheet.write(i,j+len(pos[i]), feature[i][j])
        sheet.write(i,len(pos[i])+len(feature[i]),label[i][0])
        sheet.write(i, len(pos[i]) + len(feature[i]+1), gnd[i][0])
    excel.save(r'/mnt/hgfs/transposon/myfeature1.xls')

