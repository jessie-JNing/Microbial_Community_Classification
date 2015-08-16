#!/usr/bin/env python
# encoding: utf-8

"""
Explore the properties of the feature space
1) basic information: classes, samples size, features number
2) data sparsity
3) feature appearance distribution

@author: Jessie
@email: jessie.JNing@gmail.com
"""

import pandas as pd
import math

class Feature_Exploration(object):

    def __init__(self, df):
        """
        df /DataFrame/: input features in DataFrame format
        """
        self.df = df
        self.sample_size = len(self.df.index)
        self.feature_num = len(self.df.columns)-2
        classes = self.df.loc[self.df.index, "Site"]
        self.classes = list(set(classes))

    def get_sample_size(self):

        return self.sample_size

    def get_feature_num(self):

       return self.feature_num

    def get_classes(self):

        return self.classes

    def data_sparsity(self):
        self.zero, self.non_zero = 0,0
        for col in self.df.columns[2:]:
            otu = self.df[col]
            for ele in otu:
                if ele>0:
                    self.non_zero += 1
                else:
                    self.zero += 1

    def get_zero(self):

        return self.zero

    def get_non_zero(self):

        return self.non_zero


    def otu_distribution(self):
        top = math.log(self.sample_size,2)
        num_bin = dict([(x, 0) for x in range(int(top)+2)])
        for col in self.df.columns[2:]:
            non_zero = 0
            otu = self.df[col]
            for ele in otu:
                if ele>0:
                    non_zero += 1
            log_bin = math.ceil(math.log(non_zero,2))
            num_bin[log_bin] +=1

        return num_bin


if __name__=="__main__":
    address = "/Users/Jessie/Dropbox/Pipline/OTU_Data/Teeth/otu/"
    feat_df = pd.read_csv(address + "otu_abundance.csv")
    feature_explor = Feature_Exploration(feat_df)
    print feature_explor.otu_distribution()
    

