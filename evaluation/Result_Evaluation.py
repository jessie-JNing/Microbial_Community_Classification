#!/usr/bin/env python
# encoding: utf-8

"""
Evaluate the performance of the model, including:
accuracy;
precison recall and F-score
confusion matrix

@author: Jessie
@email: jessie.JNing@gmail.com
"""

import pandas as pd
import math

class Result_Evaluation(object):

    def __init__(self, result_df):
        self.df = result_df

    def get_four_measure(self, sub_class):
        """
        sub_class /string/: target evaluated class
        """
        label = self.df.loc[self.df.index, "label"]
        prediction = self.df.loc[self.df.index, "prediction"]
        self.TP, self.FP, self.TN, self.FN = 0.0, 0.0, 0.0, 0.0
        for ele in zip(label, prediction):
            if ele[0]==sub_class and ele[1]==sub_class:
                self.TP+=1
            elif ele[0] is not sub_class and ele[1]==sub_class:
                self.FP+=1
            elif not (ele[0]==sub_class or ele[1]==sub_class):
                self.TN+=1
            elif ele[0]==sub_class and ele[1] is not sub_class:
                self.FN+=1

    def get_tp(self):

        return self.TP

    def get_fp(self):

        return self.FP

    def get_tn(self):

        return self.TN

    def get_fn(self):

        return self.FN

    def get_accuracy(self):
        label = self.df.loc[self.df.index, "label"]
        prediction = self.df.loc[self.df.index, "prediction"]
        num=0.0
        for ele in zip(label, prediction):
            if ele[0]==ele[1]:
                num+=1

        return num/len(label)


    def get_precision(self):

        return self.get_tp()/(self.get_tp() + self.get_fp())

    def get_recall(self):

        return self.get_tp()/(self.get_tp() + self.get_fn())

    def get_specificity(self):

        return self.get_tn()/(self.get_tn() + self.get_fp())

    def get_f_score(self):

        return 2*self.get_tp()/(2*self.get_tp() + self.get_fp() + self.get_fn())

    def get_mcc(self):
        upper = self.get_tp()*self.get_tn() - self.get_fn()*self.get_fp()
        lower = math.sqrt((self.get_tp() + self.get_fp())*(self.get_tp() + self.get_fn())*(self.get_tn() + self.get_fp())*(self.get_tn() + self.get_fn()))

        return upper/lower

    def get_confusion_matrix(self):

        label = self.df.loc[self.df.index, "label"]
        prediction = self.df.loc[self.df.index, "prediction"]
        classes = list(set(self.df.loc[self.df.index, "label"]))
        classes_dic = dict([(classes[i], i) for i in range(len(classes))])

        confusion_matrix = [["class"]+classes]
        for c in classes:
            actual= [c] + [0.0]*len(classes)
            confusion_matrix.append(actual)

        for ele in zip(label, prediction):
            row, col = classes_dic[ele[0]], classes_dic[ele[1]]
            confusion_matrix[row+1][col+1] += 1

        return confusion_matrix


