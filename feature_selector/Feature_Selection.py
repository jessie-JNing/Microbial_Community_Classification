#!/usr/bin/env python
# encoding: utf-8

"""
Three feature selection methods
1) information gain
2) Chi-square

@author: Jessie
@email: jessie.JNing@gmail.com
"""

import pandas as pd
import numpy as np
from sklearn.feature_selection import SelectKBest
from sklearn.feature_selection import chi2
from sklearn.ensemble import RandomForestClassifier

class Feature_Selection(object):

    def __init__(self):
        pass

    def information_gain(self, df):
        """
        df /DataFrame/: input features
        """
        X = np.array(df.loc[df.index, df.columns[2:]])
        y = np.array(df.loc[df.index, "Site"])
        def _calIg():
            entropy_x_set = 0
            entropy_x_not_set = 0
            for c in classCnt:
                probs = classCnt[c] / float(featureTot)
                entropy_x_set = entropy_x_set - probs * np.log(probs)
                probs = (classTotCnt[c] - classCnt[c]) / float(tot - featureTot)
                entropy_x_not_set = entropy_x_not_set - probs * np.log(probs)
            for c in classTotCnt:
                if c not in classCnt:
                    probs = classTotCnt[c] / float(tot - featureTot)
                    entropy_x_not_set = entropy_x_not_set - probs * np.log(probs)
            return entropy_before - ((featureTot / float(tot)) * entropy_x_set
                                 +  ((tot - featureTot) / float(tot)) * entropy_x_not_set)

        tot = X.shape[0]
        classTotCnt = {}
        entropy_before = 0
        for i in y:
            if i not in classTotCnt:
                classTotCnt[i] = 1
            else:
                classTotCnt[i] = classTotCnt[i] + 1
        for c in classTotCnt:
            probs = classTotCnt[c] / float(tot)
            entropy_before = entropy_before - probs * np.log(probs)

        nz = X.T.nonzero()
        pre = 0
        classCnt = {}
        featureTot = 0
        information_gain = []
        for i in range(0, len(nz[0])):
            if (i != 0 and nz[0][i] != pre):
                for notappear in range(pre+1, nz[0][i]):
                    information_gain.append(0)
                ig = _calIg()
                information_gain.append(ig)
                pre = nz[0][i]
                classCnt = {}
                featureTot = 0
            featureTot = featureTot + 1
            yclass = y[nz[1][i]]
            if yclass not in classCnt:
                classCnt[yclass] = 1
            else:
                classCnt[yclass] = classCnt[yclass] + 1
        ig = _calIg()
        information_gain.append(ig)

        ig=[]
        for x in zip(df.columns[2:], np.asarray(information_gain)):
            ig.append((x[0], x[1]))
        ig_sort = sorted(ig, key=lambda x: x[1], reverse=True)

        return [x[0] for x in ig_sort]



    def chi_square(self, df):
        """
        df /DataFrame/: input features
        """
        X = np.array(df.loc[df.index, df.columns[2:]])
        y = np.array(df.loc[df.index, "Site"])
        ch2 = SelectKBest(chi2, k=3).fit(X, y)
        cs = []
        for x in zip(df.columns[2:], np.asarray(ch2.scores_)):
            cs.append((x[0], x[1]))
        cs_sort = sorted(cs, key=lambda x: x[1], reverse=True)

        return [x[0] for x in cs_sort]


    def rf_feature_permutation(self, df):
        """
        df /DataFrame/: input features
        """
        X = np.array(df.loc[df.index, df.columns[2:]])
        y = np.array(df.loc[df.index, "Site"])
        rfc = RandomForestClassifier(n_estimators=100)
        rfc.fit(X, y)
        fp = []
        for x in zip(df.columns[2:], rfc.feature_importances_):
            fp.append((x[0], x[1]))
        fp_sort = sorted(fp, key=lambda x: x[1], reverse=True)

        return [x[0] for x in fp_sort]



if __name__=="__main__":
    address = "/Users/Jessie/Dropbox/Pipline/OTU_Data/Teeth/otu/"
    otu_df = pd.read_csv(address + "Chi_Square_abundance.csv")
    feature_selector = Feature_Selection()
    print feature_selector.chi_square(otu_df)