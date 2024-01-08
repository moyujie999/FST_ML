import utils
import pandas as pd
import shap
import numpy as np
from sklearn.preprocessing import LabelEncoder
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor
from sklearn.ensemble import RandomForestClassifier
from sklearn import  svm
from math import sqrt
from sklearn.metrics import mean_absolute_error
from sklearn.metrics import mean_squared_error
from sklearn.metrics import r2_score ,explained_variance_score, accuracy_score,precision_score,roc_auc_score
import matplotlib.pyplot as plt
from sklearn.model_selection import cross_val_score
from sklearn.feature_selection import RFE
from sklearn.svm import SVR
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.svm import SVC
seed = 0
np.random.seed(seed)
from sklearn.model_selection import cross_val_score
import seaborn as sns
from scipy import stats
from sklearn.model_selection import KFold
from utils import load_ra

design_orin = pd.read_csv("/home/moyujie/lab102/Microbiome_V2/data/design.csv", index_col=0)
design_orin = design_orin.drop('FST_78P', axis=0)
design_orin = design_orin.reset_index(drop=True)
AUC_score = np.zeros((4,100))

for dep in range(0,4):
    print(dep)
    if dep == 3:
        design,ra, taxa = utils.load_data_depth(time = 11, domain= 'fungi', depth=dep , random=True)
    else:
        design,ra, taxa = utils.load_data_depth(time = 12, domain= 'fungi', depth=dep , random=True)

    X = ra.reset_index()
    dep_list = ["0-10 cm", "10-20 cm", "20-30 cm", "30-60 cm"]
    chosed = design_orin[(design_orin.Sample_Depth == dep_list[dep])].index.tolist()
    design_orin_chosed = design_orin.iloc[chosed, :]

    y = design_orin_chosed.reset_index().Fertility_Source
    for i in range(10):
        rf = RandomForestClassifier()
        kfold = KFold(n_splits=10, shuffle=True, random_state = 0)
        index = kfold.split(X)
        y_test_all = []
        y_predict_all = []

        for train_index, test_index in index:
            rf2 = RandomForestClassifier()
            X_train = X.iloc[train_index, :]
            X_test = X.iloc[test_index, :]
            y_train = y[train_index]
            y_test = y[test_index]
            # X_train, y_train, X_test, y_test = kw_test(X_train, y_train, X_test, y_test, feature_num=100)
            rf2.fit(X_train, y_train)
            y_predict = rf2.predict(X_test)
            y_test_all.extend(y_test)
            y_predict_all.extend(y_predict)
        # print(accuracy_score(y_test_all,y_predict_all))
        AUC_score[dep,i]= roc_auc_score(pd.get_dummies(y_test_all), pd.get_dummies(y_predict_all),average='macro' ,multi_class= 'ovo')
        print(i)


