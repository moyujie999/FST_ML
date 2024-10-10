import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import KFold
from sklearn.metrics import roc_auc_score
from utils import load_ra

# Set seed for reproducibility
seed = 0
np.random.seed(seed)

design_orin = pd.read_csv("/home/moyujie/lab102/Microbiome_V2/data/design.csv", index_col=0)
design_orin = design_orin.drop('FST_78P', axis=0)
design_orin = design_orin.reset_index(drop=True)
AUC_score = np.zeros((4,10))

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
        kfold = KFold(n_splits=10, shuffle=True, random_state = 0)
        index = kfold.split(X)
        y_test_all = []
        y_predict_all = []

        for train_index, test_index in index:
            rf = RandomForestClassifier()
            X_train = X.iloc[train_index, :]
            X_test = X.iloc[test_index, :]
            y_train = y[train_index]
            y_test = y[test_index]
            # X_train, y_train, X_test, y_test = kw_test(X_train, y_train, X_test, y_test, feature_num=100)
            rf.fit(X_train, y_train)
            y_predict = rf.predict(X_test)
            y_test_all.extend(y_test)
            y_predict_all.extend(y_predict)
        # print(accuracy_score(y_test_all,y_predict_all))
        AUC_score[dep,i]= roc_auc_score(pd.get_dummies(y_test_all), pd.get_dummies(y_predict_all),average='macro' ,multi_class= 'ovo')
        print(i)


