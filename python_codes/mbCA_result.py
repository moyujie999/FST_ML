import utils
import numpy as np
from sklearn.metrics.pairwise import pairwise_distances
from sklearn.preprocessing import MinMaxScaler
import pandas as pd
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import KFold
from sklearn.metrics import r2_score,explained_variance_score
import shap
from sklearn.decomposition import PCA
import warnings
warnings.filterwarnings("ignore")
from sklearn.manifold import TSNE

import matplotlib.pyplot as plt



np.random.seed(0)


def DLPEstimation(X,Y):
    X_order = np.empty_like(X)
    Y_order = np.empty_like(Y)
    for i in range(X.shape[0]):
        X_order[i, :] = np.argsort(X[i, :])
        Y_order[i, :] = np.argsort(Y[i, :])
    dlp_all = 0
    start = int(X.shape[0]*0.05)
    end = int(X.shape[0]*0.5)

    for i in range(X.shape[0]):
        for j in range(start, end):
            intersect = len(np.intersect1d(X_order[i, 0:j], Y_order[i, 0:j]))
            dlp_all += intersect / j
    dlp = dlp_all / ((end-start) * X.shape[0])
    return dlp

design_org = pd.read_csv("/home/moyujie/lab102/Microbiome_V2/data/design.csv", index_col=0)
design_org = design_org.drop('FST_78P', axis=0)
design_org = design_org.reset_index(drop=True)
ra,taxa = utils.load_ra(domain= 'fungi', time = 0)

design = utils.load_design(drop= False, one_hot=False)
design['random'] = np.random.rand(design.shape[0]).tolist()
d_matrix = pairwise_distances(ra.values, metric="braycurtis")

# perplexity = 83 for prokaryotes and 90 for fungi

tsne = TSNE(metric="braycurtis", n_components=2, perplexity=90,random_state=0)
# tsne = TSNE(n_components=2)
ra_tsne = tsne.fit_transform(ra)

# plt.scatter(ra_tsne[:,0], ra_tsne[:,1])
tsne_distance = pairwise_distances(ra_tsne)
indices = np.triu_indices(tsne_distance.shape[0], k=1)
# #  Shepard diagram
# flattened_umap = umap_distance[indices]
scaler = MinMaxScaler(feature_range=(-1, 1))
normalized_ra_tsne = scaler.fit_transform(ra_tsne)
normalized_tsne_distance = pairwise_distances(normalized_ra_tsne)

flattened_tsne = tsne_distance[indices]
# flattened_pcoa = pcoa_distance[indices]
flattened_d = d_matrix[indices]
flattened_tsne_normed = normalized_tsne_distance[indices]

# plt.scatter(flattened_d, flattened_tsne)

# correlation, p_value = pearsonr(flattened_d, flattened_tsne)
# print(correlation)

explained_variance_ordination = utils.explained_variance_z_score(flattened_d, flattened_tsne)
dlp_est = DLPEstimation(d_matrix, tsne_distance)

pca = PCA()
pca.fit(ra_tsne)
principal_components = pca.components_
explained_variance_ratio = pca.explained_variance_ratio_
projected_data = pca.transform(ra_tsne)


ordination = pd.DataFrame(projected_data, columns=['x','y'])
data = pd.concat([ordination, design_org], axis=1)

# ############ Fertility Source#####################
# x = data['x']
# y = data['y']
# categories = data['Fertility_Source']
# color_map = {'Synthetic Fertilizer': 'red', 'Manure': 'blue', 'Legume': 'green'}
# for cat, color in color_map.items():
#     plt.scatter(x[categories == cat], y[categories == cat], color=color, label=cat)

########### Depth#####################
x = data['x']
y = data['y']
categories = data['Sample_Depth']
color_map = {'0-10 cm': 'red', '10-20 cm': 'blue', '20-30 cm': 'green', '30-60 cm': 'gray'}
for cat, color in color_map.items():
    plt.scatter(x[categories == cat], y[categories == cat], color=color, label=cat)


def ordination_norm(matrix):
    min_val = np.min(matrix)
    max_val = np.max(matrix)
    normalized_matrix = (matrix - min_val) / (max_val - min_val)
    normalized_matrix = 2 * normalized_matrix - 1
    return normalized_matrix


ra_tsne = ordination_norm(ra_tsne)
rf = RandomForestRegressor(max_features = 'sqrt')
y_0 = ra_tsne[:,0]
y_0 = np.array(y_0).reshape(-1)
y_1 = ra_tsne[:,1]
y_1 = np.array(y_1).reshape(-1)
X = design


kfold = KFold(n_splits=10, shuffle=True)
index = kfold.split(X)
y_test_all_0 = []
y_predict_all_0 = []
shapley_value_0 = np.empty((0, 12))
y_test_all_1 = []
y_predict_all_1 = []
shapley_value_1 = np.empty((0, 12))

for train_index, test_index in index:


    rf.fit(X.iloc[train_index, :], y_0[train_index])
    y_predict_0 = rf.predict(X.iloc[test_index, :])
    explainer = shap.TreeExplainer(rf)
    shap_values_ = explainer.shap_values(X.iloc[test_index, :])
    shap_values_ = np.abs(shap_values_)
    # shap_values_ = np.sum(shap_values_, axis=0)
    shapley_value_0 = np.vstack((shapley_value_0, shap_values_))
    y_test_all_0.extend(y_0[test_index])
    y_predict_all_0.extend(y_predict_0)

    rf.fit(X.iloc[train_index, :], y_1[train_index])
    y_predict_1 = rf.predict(X.iloc[test_index, :])
    explainer = shap.TreeExplainer(rf)
    shap_values_ = explainer.shap_values(X.iloc[test_index, :])
    shap_values_ = np.abs(shap_values_)
    # shap_values_ = np.sum(shap_values_, axis=0)
    shapley_value_1 = np.vstack((shapley_value_1, shap_values_))
    y_test_all_1.extend(y_1[test_index])
    y_predict_all_1.extend(y_predict_1)

score_0 = explained_variance_score(y_test_all_0, y_predict_all_0)
score_1 = explained_variance_score(y_test_all_1, y_predict_all_1)


shapley_value = np.sqrt(shapley_value_1**2 + shapley_value_0**2)
shapley_value_sum = np.sum(shapley_value, axis=0)

sorted_indices = np.argsort(shapley_value_sum)
# shapley_value_sum_0 = np.sum(np.array(shapley_value_0), axis=0)
# shapley_value_sum_1 = np.sum(np.array(shapley_value_1), axis=0)
sorted_columns = np.take(design.columns, sorted_indices)


explained_variance_all = explained_variance_ordination * (explained_variance_ratio[0]*score_0 + explained_variance_ratio[1]*score_1)

print(explained_variance_ordination *explained_variance_ratio[0]*score_0)
print(explained_variance_ordination *explained_variance_ratio[1]*score_1)

plt.legend()
plt.title("Scatter Plot exlained variance {:.3f}".format(explained_variance_all))
plt.xlabel("First t-SNE Axis {:.3f}".format(explained_variance_ratio[0]))
plt.ylabel("Second t-SNE Axis {:.3f}".format(explained_variance_ratio[1]))
