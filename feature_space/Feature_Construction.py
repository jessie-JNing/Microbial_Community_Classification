#!/usr/bin/env python
# encoding: utf-8

"""
Create four different feature space:
1) OTU
2) Clade
3) Function
4) Hybrid

@author: Jessie
@email: jessie.JNing@gmail.com
"""

import glob
import pandas as pd
import dendropy

from util.Sample_Mapping import Sample_Mapping


class Feature_Construction(object):

    def __init__(self):
        pass

    def create_otu_features(self, labels, otu_directory):
        """
        labels /list/: samples with the labels were used to classify
        otu_directory /string/: the directory contains picked otu files of each samples
        """

        # collect sample information
        collect_samples = {}
        appear_otus = {"Site":0, "Taxon":1}
        appear_otus_list = ["Site", "Taxon"]
        index = 2
        for label in labels:
            for samples in glob.glob(otu_directory + label + "/*_pick_otu/*_otus.txt"):
                sample_otu = {}
                sample_id = '+'.join([samples.split("/")[-1].split("_")[0], label])
                # read OTU file
                with open(samples) as file:
                    while True:
                        otu = file.readline().split()
                        if len(otu)<2:
                            break
                        else:
                            sample_otu[otu[0]] = len(otu)-1
                            if not appear_otus.has_key(otu[0]):
                                appear_otus[otu[0]]=index
                                appear_otus_list.append(otu[0])
                                index += 1
                collect_samples[sample_id] = sample_otu

        # write features into a dataframe
        otu_df = pd.DataFrame(columns=appear_otus_list)
        row_index = 0
        for sample in collect_samples:
            sample_row = sample.split("+") + [0 for x in range(len(appear_otus)-2)]
            sample_dic = collect_samples[sample]
            for item in sample_dic:
                sample_row[appear_otus[item]] = sample_dic[item]
            otu_df.loc[row_index]=sample_row
            row_index += 1

        return otu_df

    def create_clade_features(self, clade_df, tree_file, tree_format="newick"):
        """
        clade_df /DataFrame/: OTU feature in DataFrame formate
        tree /string/: directory of the phylogenetic tree
        tree_format /string/: specify the tree file format, eg. newick, nexus and etc
        """
        # read the tree file, if it's unrooted, choose the midpoint as root
        tree = dendropy.Tree()
        tree.read_from_path(tree_file, tree_format)

        # navigate the tree and calculate the clade features
        clade_counter = 0
        for node in tree.postorder_node_iter():
            if not node.is_leaf():
                node.label= "U_" + str(clade_counter)
                node_list = []
                for children in node.child_nodes():
                    node_list.append(clade_df[children.get_node_str()])
                node_sum = [sum(x) for x in zip(node_list[0], node_list[1])]
                clade_df[node.get_node_str()] = node_sum
                clade_counter += 1

        return clade_df

    def create_function_features(self, function_file):
        """
        function_file /string/: directory of PICRUSt generated file
        """
        sample_mapper = Sample_Mapping()
        mapping_dic = sample_mapper.nine_way_mapping()

        function_df = pd.DataFrame()
        with open(function_file) as f:
            f.readline()
            sample_id = f.readline().strip('\n').strip('\r').split('\t')[1:]
            print sample_id
            sample_lable = [mapping_dic[x] for x in sample_id if len(x)>5]
            function_df["Site"] = sample_lable
            function_df["Taxon"] = sample_id[:len(sample_lable)]

            while True:
                ko = f.readline().strip('\n').strip('\r').split('\t')
                if len(ko)<2:
                    break
                else:
                    ko_id = ko[0]
                    ko_value = [float(x) for x in ko[1:] if len(x)>0]
                    if sum(ko_value)>0:
                        function_df[ko_id] = ko_value

        return function_df

    def create_hybrid_features(self, clade_df, clade_feat, function_df, fun_feat):
        """
        clade_df /DataFrame/: clade feature in DataFrame format
        clade_feat /list/: clade features added to hybrid
        function_df /DataFrame/: function feature in DataFrame format
        fun_feat /list/: function features added to hybrid
        """
        header_df = clade_df.loc[clade_df.index, ["Site", "Taxon"]]
        sub_clade_df = clade_df.loc[clade_df.index, clade_feat]
        sub_fun_df = function_df.loc[function_df.index, fun_feat]

        hybrid_df = pd.concat([header_df, sub_clade_df, sub_fun_df], axis=1)

        return hybrid_df


