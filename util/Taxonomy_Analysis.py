#!/usr/bin/env python
# encoding: utf-8

"""
Analyze the proportion of of taxa at different taxonomy level.

@author: Jessie
@email: jessie.JNing@gmail.com
"""
import pandas as pd

from util import Sample_Mapping


class Taxonomy_Analysis(object):

    def __init__(self):
        pass

    def calculate_proportion(self, input_file):
        """
        input_file /string/: address of the input taxonomy level file
        """
        self.proportion_dic = {}
        label_index = {}
        sample_mapper = Sample_Mapping()
        mapping = sample_mapper.nine_way_mapping()

        with open(input_file) as f:
            # initialize the header file
            header = f.readline().strip('\n').strip('\r').split()[1:]
            for i in xrange(len(header)):
                label = mapping[header[i]]
                if label_index.has_key(label):
                    label_index[label].append(i)
                else:
                    label_index[label] = [i]

            # calculate the proportion
            while True:
                taxon = f.readline().split()
                if len(taxon)<2:
                    break
                else:
                    taxon_name = taxon[0]
                    taxon_value = taxon[1:]
                    taxa_label = {}
                    for item in label_index:
                        item_list = [float(taxon_value[x]) for x in label_index[item]]
                        taxa_label[item] = item_list
                    self.proportion_dic[taxon_name] = taxa_label

    def get_proportion(self):

        return self.proportion_dic

    def sum_taxa(self, abundance_thred, labels=None):
        """
        abundance_thred /float/: taxa that larger than the abundance_thred to be shown
        labels /list/: labels to be calculated
        """
        propor_df = pd.DataFrame()

        for taxon in self.proportion_dic:
            taxon_prop = self.proportion_dic[taxon]
            taxon_list = []
            if labels:
                for label in labels:
                    taxon_mean= sum(taxon_prop[label])/float(len(taxon_prop[label]))
                    taxon_list.append(taxon_mean)
            else:
                for label in taxon_prop.keys():
                    taxon_mean= sum(taxon_prop[label])/float(len(taxon_prop[label]))
                    taxon_list.append(taxon_mean)
                sum_prop = sum(taxon_list)
                if sum_prop>=abundance_thred:
                    propor_df[taxon] = taxon_list
        if labels:
            propor_df.index = labels
        else:
            propor_df.index = taxon_prop.keys()

        return propor_df


if __name__=="__main__":
    address = "/Users/Jessie/Dropbox/Pipline/OTU_Data/Teeth/tax_summarize/"
    TA = Taxonomy_Analysis()
    TA.calculate_proportion(address + "otu_table_tax_L2.txt")
    TA.sum_taxa(0.001)