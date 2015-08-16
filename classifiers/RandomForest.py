#!/usr/bin/env python
# encoding: utf-8

"""
Random Forest classifier

@author: Jessie
@email: jessie.JNing@gmail.com
"""


import pandas as pd
from sklearn.svm import SVC
import numpy as np
from sklearn.cross_validation import KFold
from sklearn.ensemble import RandomForestClassifier

class RandomForest(object):

    def __init__(self):
        pass


    def cross_validation(self, dataset, out_put, mtry=500, fold_k=5, max_feat="auto"):
        """
        df /DataFrame/: input feature space
        out_put /string/: directory saving the cross-valiation result
        """
        # construct input data
        attribute = np.array(dataset.loc[dataset.index, dataset.columns[2:]])
        label = np.array(dataset.loc[dataset.index, "Site"])
        sample_id = np.array(dataset.loc[dataset.index, "Taxon"])

        # do k-fold cross-validation

        clf = RandomForestClassifier(n_estimators=mtry, bootstrap=False, max_features=max_feat, random_state=0,criterion='gini')
        cv_result = {}
        for train, test in KFold(n=len(label),n_folds=fold_k, shuffle=False):
            clf.fit(attribute[train], label[train])
            predict = clf.predict(attribute[test])
            for x in zip(sample_id[test], label[test], predict):
                cv_result[x[0]] = [x[1], x[2]]

        # write result to dataframe and files as original order
        cv_result_df = pd.DataFrame(columns=["sample_id", "label", "prediction"])
        for id in sample_id:
            result = [id]
            result.extend(cv_result[id])
            cv_result_df.loc[id] = result
        cv_result_df.to_csv(out_put, index=False)


    def train_process(self, trainset, mtry=500, max_feat="auto"):
        """
        trainset /DataFrame/: train feature space
        """
        # construct input data
        train_attribute = np.array(trainset.loc[trainset.index, trainset.columns[2:]])
        train_label = np.array(trainset.loc[trainset.index, "Site"])

        # train random forest model
        clf = RandomForestClassifier(n_estimators=mtry, bootstrap=False, max_features=max_feat, random_state=0,criterion='gini')
        clf.fit(train_attribute, train_label)

        return clf


    def test_process(self, estimator, test_set):
        """
        estimator /object/: grid search result object
        outtest_set /DataFrame/: test feature space
        """
        # construct input data
        test_attribute = np.array(test_set.loc[test_set.index, test_set.columns[2:]])
        test_label = np.array(test_set.loc[test_set.index, "Site"])
        test_ample_id = np.array(test_set.loc[test_set.index, "Taxon"])

        # do prediction and save results into a DataFrame
        cv_result_df = pd.DataFrame(columns=["sample_id", "label", "prediction"])
        clf= estimator.best_estimator_
        predict = clf.predict(test_attribute)
        for x in zip(test_ample_id, test_label, predict):
            result = [x[0],x[1], x[2]]
            cv_result_df.loc[x[0]] = result

        return cv_result_df



if __name__=="__main__":
    address = "/Users/Jessie/Dropbox/Pipline/OTU_Data/Teeth/otu/"
    otu_df = pd.read_csv(address + "Chi_Square_abundance.csv")
    rf_classifier = RandomForest()
    rf_classifier.cross_validation(otu_df, address + "testrf.csv", fold_k=5)
