#!/usr/bin/env python
# encoding: utf-8

"""
Ensemble method with SVM and random forest as base classifiers

@author: Jessie
@email: jessie.JNing@gmail.com
"""


import pandas as pd
import numpy as np
from sklearn.cross_validation import KFold
from classifiers.SVM import svm
from classifiers.RandomForest import RandomForest


class Ensemble(object):

    def __init__(self):
        pass


    def create_training_prob(self, train_set, test_set, mtry=500, fold_k=5, max_feat="auto"):
        """
        train_set /DataFrame/: data set for training
        test_set /DataFrame/: data set for testing
        mtry /int/: number of tree in the random forest
        fold_k /int/: number of folds in cross-validation
        max_featu /int/: number of features in each tree in the random forest
        """
        # initialize the training and testing set
        train_labels, test_labels = np.array(train_set.loc[train_set.index, "Site"]), np.array(test_set.loc[test_set.index, "Site"])
        label_set = list(set(train_labels))
        train_sample_id, test_sample_id = np.array(train_set.loc[train_set.index, "Taxon"]),np.array(test_set.loc[test_set.index, "Taxon"])

        train_gen = pd.DataFrame(columns=["Site", "Taxon"] + ["svm+"+x for x in label_set] + ["rf+"+x for x in label_set])
        test_gen = pd.DataFrame(columns=["Site", "Taxon"] + ["svm+"+x for x in label_set] + ["rf+"+x for x in label_set])

        # SVM section
        svm_classifier = svm()
        svm_best_estimator = svm_classifier.get_best_estimator(train_set)
        svm_result_train, svm_result_test = self.svm_repeat(svm_classifier, svm_best_estimator, test_set)

        # RF section
        rf_classifier = RandomForest()
        rf_result_train = self.rf_repeat(rf_classifier, train_set, test_set, repeat_no=100)

        for id in train_sample_id:
            svm_rt = svm_result_train[id]
            rf_rt = rf_result_train[id]
            svm_list = [svm_rt[x] if svm_rt.has_key(x) else 0 for x in label_set]
            rf_list = [rf_rt[x] if rf_rt.has_key(x) else 0 for x in label_set]
            df_row = [svm_rt["label"], id] + svm_list + rf_list
            train_gen.loc[id] = df_row

        for id in test_sample_id:
            svm_rt = svm_result_train[id]
            rf_rt = rf_result_train[id]
            svm_list = [svm_rt[x] if svm_rt.has_key(x) else 0 for x in label_set]
            rf_list = [rf_rt[x] if rf_rt.has_key(x) else 0 for x in label_set]
            df_row = [svm_rt["label"], id] + svm_list + rf_list
            test_gen.loc[id] = df_row

        return train_gen, test_gen

    # methods of repeating SVM classification -- base classifier
    def svm_repeat(svm_classifier_obj, trained_estimator, test_s, repeat_no=100):
        repeat_train_result, repeat_test_result = {}, {}
        for i in range(repeat_no):
            each_train_result = svm_classifier_obj.cross_validation(trained_estimator, shuffle_op=True)
            each_test_result = svm_classifier_obj.test_process(trained_estimator, test_s)

            for row in each_train_result.index:
                sample_result = list(each_train_result.loc[row])
                sample_id = sample_result[0]
                sample_label = sample_result[1]
                sample_prediction = sample_result[-1]
                if repeat_train_result.has_key(sample_id):
                    if repeat_train_result[sample_id].has_key(sample_prediction):
                        repeat_train_result[sample_id][sample_prediction] +=1
                    else:
                        repeat_train_result[sample_id][sample_prediction] = 1
                else:
                    repeat_train_result[sample_id] = {"label":sample_label,sample_prediction:1}

            for row in each_test_result.index:
                sample_result = list(each_test_result.loc[row])
                sample_id = sample_result[0]
                sample_label = sample_result[1]
                sample_prediction = sample_result[-1]
                if repeat_test_result.has_key(sample_id):
                    if repeat_test_result[sample_id].has_key(sample_prediction):
                        repeat_test_result[sample_id][sample_prediction] +=1
                    else:
                        repeat_test_result[sample_id][sample_prediction] = 1
                else:
                    repeat_test_result[sample_id] = {"label":sample_label,sample_prediction:1}

        return repeat_train_result, repeat_test_result

    # methods of repeating RF classification -- base classifier
    def rf_repeat(rf_classifier_obj,train_s, test_s, repeat_no=100):
        repeat_train_result, repeat_test_result = {}, {}
        for i in range(repeat_no):
            each_train_result = rf_classifier_obj.cross_validation(train_s, shuffle_op=True)
            each_test_result = rf_classifier_obj.test_process(rf_classifier_obj, test_s)

            for row in each_train_result.index:
                sample_result = list(each_train_result.loc[row])
                sample_id = sample_result[0]
                sample_label = sample_result[1]
                sample_prediction = sample_result[-1]
                if repeat_train_result.has_key(sample_id):
                    if repeat_train_result[sample_id].has_key(sample_prediction):
                        repeat_train_result[sample_id][sample_prediction] +=1
                    else:
                        repeat_train_result[sample_id][sample_prediction] = 1
                else:
                    repeat_train_result[sample_id] = {"label":sample_label,sample_prediction:1}

            for row in each_test_result.index:
                sample_result = list(each_test_result.loc[row])
                sample_id = sample_result[0]
                sample_label = sample_result[1]
                sample_prediction = sample_result[-1]
                if repeat_test_result.has_key(sample_id):
                    if repeat_test_result[sample_id].has_key(sample_prediction):
                        repeat_test_result[sample_id][sample_prediction] +=1
                    else:
                        repeat_test_result[sample_id][sample_prediction] = 1
                else:
                    repeat_test_result[sample_id] = {"label":sample_label,sample_prediction:1}

        return repeat_train_result, repeat_test_result


    # use random forest as an ensemble classifier
    def ensemble_cross_validation(self, dataset, fold_k = 5):
        """
        dataset /DataFrame/: data set
        fold_k /int/: number of folds in cross-validation
        """

        ensemble_result = pd.DataFrame(columns=["sample_id", "label", "prediction"])
        label = np.array(dataset.loc[dataset.index, "Site"])
        for train, test in KFold(n=len(label),n_folds=fold_k, shuffle=False):
            train_set = dataset.loc[train, dataset.columns]
            test_set = dataset.loc[test, dataset.columns]

            train_att, test_att = self.create_training_prob(train_set, test_set)
            ensemble_model = RandomForest()
            ensemble_estimator = ensemble_model.train_process(train_att)
            ensemble_sub_result = ensemble_model.test_process(ensemble_estimator, test_att)
            ensemble_result = pd.concat([ensemble_result, ensemble_sub_result])

        return ensemble_result