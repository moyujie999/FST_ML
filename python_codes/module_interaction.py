import utils
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import KFold
from sklearn.metrics import r2_score
import shap
import warnings
warnings.filterwarnings("ignore")
from sklearn.preprocessing import MinMaxScaler

depth = 0


ra_b,taxa_b = utils.load_ra(domain= 'bacteria', time = 0)
ra_f,taxa_f = utils.load_ra(domain= 'fungi', time = 0)
design = utils.load_design(drop= False)
env = utils.load_design(drop= False)
env = env[env.Sample_Depth == depth]
ra = pd.concat([ra_b, ra_f], axis = 1)

chosed = design[(design.Sample_Depth != depth)].index.tolist()

ra = ra.drop(chosed, axis=0)

module_legu_0 = pd.read_csv("/home/moyujie/lab102/Microbiome_V2/data/module_data/module_legu_dep0.csv", index_col=0)
module_syn_0 = pd.read_csv("/home/moyujie/lab102/Microbiome_V2/data/module_data/module_syn_dep0.csv", index_col=0)
module_legu_1 = pd.read_csv("/home/moyujie/lab102/Microbiome_V2/data/module_data/module_legu_dep1.csv", index_col=0)
module_syn_1 = pd.read_csv("/home/moyujie/lab102/Microbiome_V2/data/module_data/module_syn_dep1.csv", index_col=0)
module_legu_2 = pd.read_csv("/home/moyujie/lab102/Microbiome_V2/data/module_data/module_legu_dep2.csv", index_col=0)
module_syn_2 = pd.read_csv("/home/moyujie/lab102/Microbiome_V2/data/module_data/module_syn_dep2.csv", index_col=0)



legu_0 = module_legu_0.index.tolist()
syn_0 = module_syn_0.index.tolist()
legu_1 = module_legu_1.index.tolist()
syn_1 = module_syn_1.index.tolist()
legu_2 = module_legu_2.index.tolist()
syn_2 = module_syn_2.index.tolist()

syn_list = [syn_0,syn_1,syn_2]
legu_list = [legu_0,legu_1,legu_2]


legu = list(set(legu_0).union(set(legu_1)).union(set(legu_2)))
syn = list(set(syn_0).union(set(syn_1)).union(set(syn_2)))


ra_legu = ra[legu]
ra_syn = ra[syn]


ra = ra[syn_list[depth]]

print(len(list(ra)))
score_all1 = []
shap_list = []



scaler = MinMaxScaler(feature_range=(0, 1))
ra = ra * 1000 # Prevent model failure due to too small data
reads = list(ra)

for i in range(len(reads)):
    print(i)
    y = ra.iloc[:, i]
    y = np.array(y).reshape(-1)
    X = ra.drop(reads[i], axis=1)

    kfold = KFold(n_splits=10, shuffle=True, random_state = 0)
    index = kfold.split(X)
    y_test_all = []
    y_predict_all = []
    shapley_value = []
    for train_index, test_index in index:
        rf = RandomForestRegressor( random_state = 0)

        rf.fit(X.iloc[train_index, :], y[train_index])
        y_predict = rf.predict(X.iloc[test_index, :])
        explainer = shap.TreeExplainer(rf)
        shap_values_ = explainer.shap_values(X.iloc[test_index, :])
        shap_values_ = np.abs(shap_values_)
        shap_values_ = np.sum(shap_values_, axis=0)
        shapley_value.append(shap_values_)
        y_test_all.extend(y[test_index])
        y_predict_all.extend(y_predict)

    score1 = r2_score(y_test_all, y_predict_all)
    score_all1.append(score1)
    shapley_value = np.sum(np.array(shapley_value), axis=0)
    shap_list.append(shapley_value)


score_all1 = np.array(score_all1)
score_mean = np.mean(score_all1)

r_score = score_all1[np.where(score_all1>0.2)]
ratio_abio = r_score.shape[0]/score_all1.shape[0]
