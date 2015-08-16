#!/usr/bin/env python
# encoding: utf-8

"""
Map labels and sample ids from the mapping file
1) basic information: classes, samples size, features number
2) data sparsity
3) otu abundance distribution
4) feature appearance distribution

@author: Jessie
@email: jessie.JNing@gmail.com
"""

class Sample_Mapping(object):

    def __init__(self):
        pass

    def nine_way_mapping(self):
        mapping_dic = {}
        with open("/Users/Jessie/Documents/Learning/CS6802/Project/Oral_Cavity/Mapping_35.txt") as f:
            for i in range(2):
                header = f.readline()
            while True:
                sample = f.readline().split()
                if len(sample)<2:
                    break
                else:
                    sample_id = sample[0]
                    sample_lable = sample[1]
                    mapping_dic[sample_id] = sample_lable
        return mapping_dic


    def four_way_mapping(self):
        mapping_dic = {}
        with open("/Users/Jessie/Documents/Learning/CS6802/Project/Oral_Cavity/Mapping_35.txt") as f:
            for i in range(2):
                header = f.readline()
            while True:
                sample = f.readline().split()
                if len(sample)<2:
                    break
                else:
                    sample_id = sample[0]
                    sample_lable = sample[-1]
                    mapping_dic[sample_id] = sample_lable
        return mapping_dic


