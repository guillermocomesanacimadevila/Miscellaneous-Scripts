#!/usr/bin/env python3
#train.py
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from joblib import dump
import shap
from sklearn.neighbors import KNeighborsClassifier # will probs be good -> BClas
from sklearn.linear_model import LogisticRegression # none / l1 / l2
from sklearn.svm import SVC # rbf kernel
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report, roc_auc_score, log_loss, confusion_matrix, precision_recall_curve
from sklearn.model_selection import train_test_split, RepeatedStratifiedKFold, GridSearchCV
#from imblearn.over_sampling import SMOTE
# deep learning stuff
#import torch
#from torch import nn

# models/${model}.py
# model -> make it CLI arg in nextflow
# KNN
# SVR-RBF
# LASSO LOGISTIC
# RIDGE LOGISTIC
# UNPENALISED LR
# XGBOOST
# MLP
# TRANSFORMER

# for utils.py
# Balancing function (SMOTE)
# Grid searches
# Stratified K-fold CV (+ CIs)
# Save models as .pkl file
# Save per-model outputs as .csv

try:
    from xgboost import XGBClassifier
    HAS_XGB = True
except:
    HAS_XGB = False

# __init__.py
# build_feature_matrix.py
# utils.py
# train.py
# predict.py
# models/LR (l1,l2,no penalty)
# models/XGBoost
# models/MLP
# models/Transformer

# only designed for binary class -> for now
def decide_smote(matrix, target_col):
    df = pd.DataFrame(matrix)
    def calc_imbalance_ratio():
        # between < 40% or greater than 60%
        counts["Counts"] = dict(df[target_col].value_counts().sort_index())

    counts = df[target_col].value_counts().sort_index()
    # if counts or <40 -> SMOTE during training
    #smote = SMOTE(random_state=0)
    # smote.fit()

    return counts

lr = LogisticRegression(
    C=1.0,
    random_state=0,
    solver='liblinear',
    penalty='l1'
)

xgb = XGBClassifier(
    n_estimators=100,
    n_jobs=-1,
    random_state=0
)

rf = RandomForestClassifier(
    n_estimators=100,
    n_jobs=-1,
    random_state=0
)

svc = SVC(
    kernel='rbf',
    gamma='auto',
    random_state=0
)

knn = KNeighborsClassifier(
    n_neighbors=5,
    n_jobs=-1,
)

# check gradients
feature_matrix = pd.read_csv("../data/ml/H37Rv_INH_presence_absence.csv")
feature_matrix = feature_matrix.rename(columns={"Unnamed: 0": "Sample"})
print(decide_smote(feature_matrix, "Label"))


# we need to ensure to split evenly across traina nd tets
X = feature_matrix.drop(["Label", "Sample"], axis=1)
y = feature_matrix["Label"]
rskf = RepeatedStratifiedKFold(n_splits=5, n_repeats=5, random_state=0)

X_train, X_test, y_train, y_test = train_test_split(
    X,
    y,
    test_size=0.2,
    random_state=0
)

gs = GridSearchCV(
    lr, {'C': [0.1, 1, 10, 100]}
) # 1 == best
gs.fit(X_train, y_train)
lr.fit(X_train, y_train)
y_pred = lr.predict(X_test)
report = classification_report(y_test, y_pred)
pr_auc = roc_auc_score(y_test, y_pred)
print("AUC: ", pr_auc)
print(report)
print(f"BCE_loss LR (l1 penalty): {log_loss(y_test, y_pred)}")

# xgb
xgb.fit(X_train, y_train)
y_pred = xgb.predict(X_test)
report = classification_report(y_test, y_pred)
pr_auc = roc_auc_score(y_test, y_pred)
print("XGB AUC: ", pr_auc)
print(report)
print(f"BCE_loss XGB: {log_loss(y_test, y_pred)}")
cm = confusion_matrix(y_test, y_pred)
#pr = precision_recall_curve(y_test, y_pred)

# rf
rf.fit(X_train, y_train)
y_pred = rf.predict(X_test)
report = classification_report(y_test, y_pred)
pr_auc = roc_auc_score(y_test, y_pred)
print("RF AUC: ", pr_auc)
print(report)
print(f"BCE_loss LR (RF): {log_loss(y_test, y_pred)}")

# rbf-svc
svc.fit(X_train, y_train)
y_pred = svc.predict(X_test)
report = classification_report(y_test, y_pred)
pr_auc = roc_auc_score(y_test, y_pred)
print("SVM AUC: ", pr_auc)
print(report)
print(f"BCE_loss LR (SVM): {log_loss(y_test, y_pred)}")

# knn
knn.fit(X_train, y_train)
y_pred = knn.predict(X_test)
report = classification_report(y_test, y_pred)
pr_auc = roc_auc_score(y_test, y_pred)
print("KNN AUC: ", pr_auc)
print(report)
print(f"BCE_loss LR (KNN): {log_loss(y_test, y_pred)}")

# RF does not perform well when X \in {0, 1}
#plt.figure()
#sns.heatmap(cm, annot=True)
#plt.show()
# print([i for i in feature_matrix.columns if i.startswith("Label")])
