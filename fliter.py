#!/usr/bin/env python
# -*- coding: utf-8 -*-

from read_function import *


def adjust(candidate_area_dict,f_human):
    """adjust start and end postion"""
    keys = sorted(candidate_area_dict.keys())
    for key in keys:
        #print key
        chr=key[0]
        start=key[1]
        end=key[2]
        while getMaxBlocks(f_human.fetch(chr, start-20, end+20))!=(start, end):
            (start, end) = getMaxBlocks(f_human.fetch(chr, start-20, end+20))
        if key!=(chr,start, end):
            candidate_area_dict[(chr,start, end)] = candidate_area_dict[key]
            del candidate_area_dict[key]

def merge(candidate_area_dict,window_size):
    keys=sorted(candidate_area_dict.keys())
    prim=()
    for key in keys:
        if candidate_area_dict[key]<4: del candidate_area_dict[key]
        if prim:
            if prim[0]==key[0] and prim[2]+window_size*0.5 >= key[1]:
                try:
                    c = (key[0], prim[1], key[2])
                    candidate_area_dict[c] = candidate_area_dict[prim] + candidate_area_dict[key]
                    #candidate_area_dict[c] = list(set(candidate_area_dict[prim] + candidate_area_dict[key]))
                    del candidate_area_dict[prim]
                    if c!=key:
                        del candidate_area_dict[key]
                        prim = c
                except KeyError:
                    continue
            else:
                prim=key
        else:
             prim = key


if __name__ == "__main__":
    candidate_list={('a',1,3):[1,2],('a',4,5):[1,3],('b',5,8):[2,3],('b',7,9):[2,4,5]}
    #print merge(candidate_list,300)