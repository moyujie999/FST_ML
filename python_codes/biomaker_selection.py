import utils
import pandas as pd
import shap
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score
import numpy as np
seed = 0
np.random.seed(seed)
from sklearn.model_selection import KFold


depth_list = [[0],[1],[2],[3],[0,1],[1,2],[2,3],[0,1,2],[1,2,3],[0,1,2,3]]
time_list = [12,12,12,12,24,24,24,36,36,48]

farming_practice = 'F'
domain = 'prokaryotes'

feature_num_list = []
feature_num_mean_score_list = []
feature_num_std_score_list = []

depth_index = 9

ra, taxa, design = utils.load_RF_task(time = time_list[depth_index],domain = domain, depth = depth_list[depth_index])
ra = ra.reset_index(drop=True)

design = design.reset_index()
X = ra
if farming_practice == 'T':
    y = design.Tillage
elif farming_practice == 'F':
    y = design.Fertility_Source
else:
    y = design.Cover_Crop

print(depth_index)
score = np.zeros((10,))

y_test_all = []
y_predict_all = []
biomarkers_frequency = {}
kfold = KFold(n_splits=10 , shuffle=True, random_state=0)
index = kfold.split(X)


for train_index, test_index in index:
    rf2 = RandomForestClassifier(n_estimators=1000, max_depth=10,max_features = 'sqrt',random_state=10)
    X_train = X.iloc[train_index, :]
    X_test = X.iloc[test_index, :]
    y_train = y[train_index]
    y_test = y[test_index]

    # X_train, y_train, X_test, y_test = kw_test(X_train, y_train, X_test, y_test, feature_num=100)
    rf2.fit(X_train, y_train)
    explainer = shap.TreeExplainer(rf2)

    shap_values = explainer.shap_values(X_train)
    shap_values_abs_sum = np.zeros_like(shap_values[0])

    for shap_values_item in shap_values:
        shap_values_abs_sum += np.abs(shap_values_item)
    shap_values_abs_sum = np.sum(shap_values_abs_sum, axis=0)
    sorted_indices = np.argsort(-shap_values_abs_sum)[0:50]


    biomraker_list = X_train.columns[sorted_indices].tolist()



    # print(len(biomraker_list))
    X_train = X_train.iloc[:, sorted_indices]
    X_test = X_test.iloc[:, sorted_indices]


    rf2 = RandomForestClassifier()
    rf2.fit(X_train, y_train)
    explainer = shap.TreeExplainer(rf2)
    shap_values = explainer.shap_values(X_train)
    shap_values_abs_sum = np.zeros_like(shap_values[0])

    for shap_values_item in shap_values:
        shap_values_abs_sum += np.abs(shap_values_item)
    shap_values_abs_sum = np.sum(shap_values_abs_sum, axis=0)

    y_predict = rf2.predict(X_test)
    y_test_all.extend(y_test)
    y_predict_all.extend(y_predict)
    j = 0
    for ASV in biomraker_list:
        if ASV in biomarkers_frequency:
            biomarkers_frequency[ASV] += shap_values_abs_sum[j]
        else:
            biomarkers_frequency[ASV] = shap_values_abs_sum[j]
        j+= 1


AUC_score = roc_auc_score(pd.get_dummies(y_test_all), pd.get_dummies(y_predict_all),average='macro' ,multi_class= 'ovo')

#

# Sort the biomarkers_frequency dictionary by its values in descending order
sorted_biomarkers = {k: v for k, v in sorted(biomarkers_frequency.items(), key=lambda item: item[1], reverse=True)}
ASV_top50 = list(sorted_biomarkers.keys())[0:50]
contribution = list(sorted_biomarkers.values())[0:50]




biomarker = taxa.loc[ASV_top50]
biomarker['Contribution'] = contribution
print(AUC_score)