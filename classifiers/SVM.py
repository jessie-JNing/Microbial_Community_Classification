#!/usr/bin/env python
# encoding: utf-8

"""
SVM classifier

@author: Jessie
@email: jessie.JNing@gmail.com
"""


import pandas as pd
from sklearn.svm import SVC
import numpy as np
from sklearn.cross_validation import KFold
from sklearn.grid_search import GridSearchCV

class svm(object):

    def __init__(self):
        pass


    def get_best_estimator(self, dataset, fold_k=5):
        """
        dataset /DataFrame/: input features
        fold_k /int/: k-fold cross-validation
        """

        # construct input data
        self.attribute = np.array(dataset.loc[dataset.index, dataset.columns[2:]])
        self.label = np.array(dataset.loc[dataset.index, "Site"])
        self.sample_id = np.array(dataset.loc[dataset.index, "Taxon"])
        self.fold_k = fold_k

        # initialize grid search parameters
        C_range = [2**i for i in range(5,15,2)]
        gamma_range = [2**i for i in range(3, -15, -2)]
        param_grid = dict(gamma=gamma_range, C=C_range)

        # use grid search to find the best parameter combination
        cv = KFold(n=len(self.label), n_folds=fold_k,shuffle=False)
        grid = GridSearchCV(SVC(probability=False, cache_size=3000), param_grid=param_grid, cv=cv,n_jobs=3)
        grid.fit(self.attribute, self.label)
        return grid

    def cross_validation(self, estimator, shuffle_op=False):
        """
        estimator /object/: grid search result object
        out_put /string/: directory saving the cross-valiation result
        """

        # do k-fold cross-validation with optimized parameters
        para = estimator.best_params_
        clf = SVC(C=para["C"], gamma=para["gamma"])
        cv_result = {}
        for train, test in KFold(n=len(self.label),n_folds=self.fold_k, shuffle=False):
            clf.fit(self.attribute[train], self.label[train])
            predict = clf.predict(self.attribute[test])
            for x in zip(self.sample_id[test], self.label[test], predict):
                cv_result[x[0]] = [x[1], x[2]]

        # write result to dataframe and files as original order
        cv_result_df = pd.DataFrame(columns=["sample_id", "label", "prediction"])
        for id in self.sample_id:
            result = [id]
            result.extend(cv_result[id])
            cv_result_df.loc[id] = result

        return cv_result_df



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


