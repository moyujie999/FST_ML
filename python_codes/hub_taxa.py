import utils
import numpy as np
import pandas as pd
import warnings
warnings.filterwarnings("ignore")

ra_b,taxa_b = utils.load_ra(domain= 'prokaryotes', time = 0)
ra_f,taxa_f = utils.load_ra(domain= 'fungi', time = 0)
# design = utils.load_design(drop= False)
design = utils.load_design(drop= False)
design = design[design.Sample_Depth == 0]
ra = pd.concat([ra_b, ra_f], axis = 1)
env = utils.load_design(drop= True)

taxa = pd.concat([taxa_b, taxa_f], axis = 0)

chosed = design[(design.Sample_Depth == 3)].index.tolist()
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

legu_list = [legu_0,legu_1,legu_2]
syn_list = [syn_0,syn_1,syn_2]

legu = list(set(legu_0).union(set(legu_1)).union(set(legu_2)))
syn = list(set(syn_0).union(set(syn_1)).union(set(syn_2)))


ra_legu = ra[legu]
ra_syn = ra[syn]

shap_legu_0 = np.load("/home/moyujie/lab102/Microbiome_V2/result_module_V4_seed/shap_legu_dep0.npy")
shap_legu_1 = np.load("/home/moyujie/lab102/Microbiome_V2/result_module_V4_seed/shap_legu_dep1.npy")
shap_legu_2 = np.load("/home/moyujie/lab102/Microbiome_V2/result_module_V4_seed/shap_legu_dep2.npy")

shap_syn_0 = np.load("/home/moyujie/lab102/Microbiome_V2/result_module_V4_seed/shap_syn_dep0.npy")
shap_syn_1 = np.load("/home/moyujie/lab102/Microbiome_V2/result_module_V4_seed/shap_syn_dep1.npy")
shap_syn_2 = np.load("/home/moyujie/lab102/Microbiome_V2/result_module_V4_seed/shap_syn_dep2.npy")

r_square_syn_0 = np.load("/home/moyujie/lab102/Microbiome_V2/result_module_V4_seed/r_square_syn_dep0.npy")
r_square_syn_1 = np.load("/home/moyujie/lab102/Microbiome_V2/result_module_V4_seed/r_square_syn_dep1.npy")
r_square_syn_2 = np.load("/home/moyujie/lab102/Microbiome_V2/result_module_V4_seed/r_square_syn_dep2.npy")

r_square_legu_0 = np.load("/home/moyujie/lab102/Microbiome_V2/result_module_V4_seed/r_square_legu_dep0.npy")
r_square_legu_1 = np.load("/home/moyujie/lab102/Microbiome_V2/result_module_V4_seed/r_square_legu_dep1.npy")
r_square_legu_2 = np.load("/home/moyujie/lab102/Microbiome_V2/result_module_V4_seed/r_square_legu_dep2.npy")



def shap_bio_preprocess(shap_list):
    ASV_num = len(shap_list[:,0])
    ASV_im = np.zeros((ASV_num,ASV_num))
    for i in range(ASV_num):
        change = list(shap_list[i, :])
        change.insert(i, 0.0)
        ASV_im[i,:] = change
    return ASV_im


def hub_taxa(r_square, shap, taxa_list,output_file):
    indices = np.where(r_square < 0.2)[0]
    shap_score_syn = shap_bio_preprocess(shap)
    shap_score_syn[:, indices] = 0
    shap_score_syn = shap_score_syn.sum(axis=1)
    order_syn = shap_score_syn.argsort()[::-1]
    shap_list_syn = shap_score_syn[order_syn]
    taxa_T = taxa.T
    taxa_syn = taxa_T[taxa_list].T
    taxa_syn_order= taxa_syn.iloc[order_syn,]
    percentage = 0.1
    num_rows_to_keep = int(len(taxa_syn_order) * percentage)
    df_top_10_percent = taxa_syn_order.head(num_rows_to_keep)
    df_top_10_percent.to_csv(output_file, index=True)
#
hub_taxa(r_square = r_square_syn_0, shap = shap_syn_0, taxa_list = syn_0, output_file="/home/moyujie/lab102/Microbiome_V2/hub_taxa/v2_seed/hub_syn_dep0.csv")
hub_taxa(r_square = r_square_syn_1, shap = shap_syn_1, taxa_list = syn_1, output_file="/home/moyujie/lab102/Microbiome_V2/hub_taxa/v2_seed/hub_syn_dep1.csv")
hub_taxa(r_square = r_square_syn_2, shap = shap_syn_2, taxa_list = syn_2, output_file="/home/moyujie/lab102/Microbiome_V2/hub_taxa/v2_seed/hub_syn_dep2.csv")

hub_taxa(r_square = r_square_legu_0, shap = shap_legu_0, taxa_list = legu_0, output_file="/home/moyujie/lab102/Microbiome_V2/hub_taxa/v2_seed/hub_legu_dep0.csv")
hub_taxa(r_square = r_square_legu_1, shap = shap_legu_1, taxa_list = legu_1, output_file="/home/moyujie/lab102/Microbiome_V2/hub_taxa/v2_seed/hub_legu_dep1.csv")
hub_taxa(r_square = r_square_legu_2, shap = shap_legu_2, taxa_list = legu_2, output_file="/home/moyujie/lab102/Microbiome_V2/hub_taxa/v2_seed/hub_legu_dep2.csv")

