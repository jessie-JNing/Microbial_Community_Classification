#!/usr/bin/env python
# encoding: utf-8

"""
Normalize features in three ways
1) OTU count
2) TSS convert
3) CSS convert

@author: Jessie
@email: jessie.JNing@gmail.com
"""
class Feature_Normalization(object):

    def __init__(self):
        pass

    def normalize_features(self, feature_df, normalization ,out_put):
        """
        feature_df /DataFrame/: all the features in a DataFrame structure
        normalization /string/: "count", "TSS", "CSS"
        out_put /string/: expected output file directory and name
        """
        if normalization=="TSS":
            print "TSS normalization"
            for i in range(len(feature_df.index)):
                sample_row = list(feature_df.loc[i])
                seq_sum = sum([float(x) for x in sample_row[2:]])
                TSS_row = sample_row[:2] + [x/seq_sum for x in sample_row[2:]]
                feature_df.loc[i] = TSS_row

        elif normalization == "CSS":
            print "CSS normalization"
            for i in range(len(feature_df.index)):
                sample_row = list(feature_df.loc[i])
                seq_sort = sorted([x for x in sample_row[2:] if x>0])
                percentile_length = int(len(seq_sort)*0.75)
                seq_sum = sum([x for x in seq_sort[:percentile_length]])
                CSS_row = sample_row[:2] + [float(x)/seq_sum for x in sample_row[2:]]
                feature_df.loc[i] = CSS_row

        feature_df.to_csv(out_put, index=False)