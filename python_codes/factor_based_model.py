import utils
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import KFold
from sklearn.metrics import r2_score
import warnings
warnings.filterwarnings("ignore")
from sklearn.preprocessing import MinMaxScaler

factor = 'bio'
domain = 'prokaryotes'
dep = 1

design,ra, taxa = utils.load_data_depth(time = 12, domain= domain, depth=dep,random=True)
reads = list(ra)
score_all1 = []
shap_list = []

scaler = MinMaxScaler(feature_range=(0, 1))
ra = ra * 1000  # Prevent model failure due to too small data
for i in range(len(reads)):
    print(i)
    y = ra.iloc[:, i]
    y = np.array(y).reshape(-1)

    if factor == 'bio':
        X = ra.drop(reads[i], axis=1)
    else:
        X = design


    kfold = KFold(n_splits=10, shuffle=True)
    index = kfold.split(X)
    y_test_all = []
    y_predict_all = []
    shapley_value = []
    for train_index, test_index in index:
        rf = RandomForestRegressor()

        rf.fit(X.iloc[train_index, :], y[train_index])
        y_predict = rf.predict(X.iloc[test_index, :])


        y_test_all.extend(y[test_index])
        y_predict_all.extend(y_predict)

    score1 = r2_score(y_test_all, y_predict_all)
    score_all1.append(score1)
    # shapley_value = np.sum(np.array(shapley_value), axis=0)
    # shap_list.append(shapley_value)


# #
def evaluate(score):
    score = np.array(score)
    score = np.abs(score)
    score = np.mean(score)
    return score
score_all1 = np.array(score_all1)
score_mean = evaluate(score_all1)

r_score = score_all1[np.where(score_all1>0.2)]
ratio_abio = r_score.shape[0]/score_all1.shape[0]