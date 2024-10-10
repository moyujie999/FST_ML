import pandas as pd
import numpy as np
from sklearn.preprocessing import LabelEncoder
from sklearn.metrics import r2_score,explained_variance_score

def load_ra(time = 0,domain = 'prokaryotes', drop = False):

    if domain == 'prokaryotes':
        ra = pd.read_csv("/home/moyujie/lab102/Microbiome_V2/data/ASVs_prokaryotes.csv", index_col= 0 )
        taxa = pd.read_csv("/home/moyujie/lab102/Microbiome_V2/data/taxa_prokaryotes.csv", index_col=0)
        ra = ra.T
        ra = ra.drop('FST_78P', axis=0)
    else:
        ra = pd.read_csv("/home/moyujie/lab102/Microbiome_V2/data/ASVs_fungi.csv", index_col=0)
        taxa = pd.read_csv("/home/moyujie/lab102/Microbiome_V2/data/taxa_fungi.csv", index_col=0)
        ra = ra.T
        ra = ra.drop('FST_78F', axis=0)


    ra = ra.reset_index(drop=True)
    ra = ra.div(ra.sum(axis=1), axis=0)
    num = (ra == 0).astype(int).sum(axis=0)
    a = np.array(np.where(np.array(num) < (96-time))).reshape(-1)
    ra = ra.iloc[:, a]
    taxa = taxa.iloc[a,:]
    return ra, taxa


def load_RF_task(time = 0,domain = 'prokaryotes', depth = [0,1,2,3]):
    if domain == 'prokaryotes':
        ra = pd.read_csv("/home/moyujie/lab102/Microbiome_V2/data/ASVs_prokaryotes.csv", index_col= 0 )
        taxa = pd.read_csv("/home/moyujie/lab102/Microbiome_V2/data/taxa_prokaryotes.csv", index_col=0)
        ra = ra.T
        ra = ra.drop('FST_78P', axis=0) # Delete the abnormal sample
    else:
        ra = pd.read_csv("/home/moyujie/lab102/Microbiome_V2/data/ASVs_fungi.csv", index_col=0)
        taxa = pd.read_csv("/home/moyujie/lab102/Microbiome_V2/data/taxa_fungi.csv", index_col=0)
        ra = ra.T
        ra = ra.drop('FST_78F', axis=0)

    design = pd.read_csv("/home/moyujie/lab102/Microbiome_V2/data/design.csv", index_col=0)
    design = design.drop('FST_78P', axis=0)
    design = design.reset_index(drop=True)
    ra = ra.reset_index(drop=True)
    ra = ra.div(ra.sum(axis=1), axis=0)

    dep_0 = design[(design.Sample_Depth == '0-10 cm')].index.tolist()
    dep_1 = design[(design.Sample_Depth == '10-20 cm')].index.tolist()
    dep_2 = design[(design.Sample_Depth == '20-30 cm')].index.tolist()
    dep_3 = design[(design.Sample_Depth == '30-60 cm')].index.tolist()

    depth_list = [dep_0,dep_1,dep_2,dep_3]
    chosed = []
    for dep in depth:
        chosed += depth_list[dep]

    ra = ra.iloc[chosed, :]
    design = design.iloc[chosed, :]
    num = (ra == 0).astype(int).sum(axis=0)
    a = np.array(np.where(np.array(num) < (ra.shape[0]-time))).reshape(-1)
    ra = ra.iloc[:, a]
    taxa = taxa.iloc[a,:]
    return ra, taxa, design



def load_design(drop = False, one_hot = True):
    '''
    # 0:Manure, 1:synthetic fertilizer, 2: legume
    # 0: reduced till, 1:till
    # 0: cover crop, no cover crop

    no one hot: 0: "Synthetic Fertilizer", 1: "Manure", 2: "Legume"

    '''
    design = pd.read_csv("/home/moyujie/lab102/Microbiome_V2/data/design.csv", index_col=0)
    design = design.drop('FST_78P', axis=0)
    design = design.reset_index(drop=True)
    if drop == True:
        design = design.iloc[:,
                 [6, 18, 21, 22, 25, 26, 27, 28, 32, 33, 34, 35, 36, 37, 38, 39, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50,
                  51, 52]]
        fertility_label = pd.get_dummies(design.Fertility_Source, prefix='Fertility')
        encoder_depth = LabelEncoder()
        encoder_management = LabelEncoder()
        encoder_tillage = LabelEncoder()
        encoder_cover = LabelEncoder()
        encoder_depth.fit(["0-10 cm", "10-20 cm", "20-30 cm", "30-60 cm"])
        encoder_tillage.fit(["Reduced Till", "Till"])
        encoder_cover.fit(["Cover Crop", "No Cover Crop"])
        design.Sample_Depth = encoder_depth.transform(design.Sample_Depth)
        design.Tillage = encoder_tillage.transform(design.Tillage)
        design.Cover_Crop = encoder_cover.transform(design.Cover_Crop)

        if one_hot == True:

            design = design.drop("Fertility_Source", axis=1)
            design = design.join(fertility_label)
        else:
            encoder_fertility = LabelEncoder()
            encoder_fertility.fit(["Synthetic Fertilizer", "Manure", "Legume"])
            design.Fertility_Source = encoder_fertility.transform(design.Fertility_Source)
        chosed = design[(design.Sample_Depth == 3)].index.tolist()
        design = design.drop(chosed, axis=0)
    else:
        design = design.iloc[:, [6, 18, 21, 22, 25, 26, 27, 33,34, 35,36]]
        fertility_label = pd.get_dummies(design.Fertility_Source, prefix = 'Fertility' )
        encoder_depth = LabelEncoder()
        encoder_management = LabelEncoder()
        encoder_tillage = LabelEncoder()
        encoder_cover = LabelEncoder()
        encoder_depth.fit(["0-10 cm", "10-20 cm", "20-30 cm", "30-60 cm"])
        encoder_tillage.fit(["Reduced Till", "Till"])
        encoder_cover.fit(["Cover Crop", "No Cover Crop"])
        design.Sample_Depth = encoder_depth.transform(design.Sample_Depth)
        design.Tillage = encoder_tillage.transform(design.Tillage)
        design.Cover_Crop = encoder_cover.transform(design.Cover_Crop)
        if one_hot == True:

            design = design.drop("Fertility_Source", axis=1)
            design = design.join(fertility_label)
        else:
            encoder_fertility = LabelEncoder()
            encoder_fertility.fit(["Synthetic Fertilizer", "Manure", "Legume"])
            design.Fertility_Source = encoder_fertility.transform(design.Fertility_Source)
    return design


def load_data_depth(time = 0,domain = 'prokaryotes', depth = 0, random = False):
    '''
    # 0:Manure, 1:synthetic fertilizer, 2: legume
    # 0: reduced till, 1:till
    # 0: cover crop, no cover crop
    '''
    design = pd.read_csv("/home/moyujie/lab102/Microbiome_V2/data/design.csv", index_col=0)
    design = design.drop('FST_78P', axis=0)
    design = design.reset_index(drop=True)
    design = design.iloc[:,
             [6, 18, 21, 22, 25, 26, 27, 28, 32, 33, 34, 35, 36, 37, 38, 39, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51,
              52]]
    fertility_label = pd.get_dummies(design.Fertility_Source, prefix='Fertility')
    encoder_depth = LabelEncoder()
    encoder_management = LabelEncoder()
    encoder_tillage = LabelEncoder()
    encoder_cover = LabelEncoder()
    encoder_depth.fit(["0-10 cm", "10-20 cm", "20-30 cm", "30-60 cm"])
    encoder_tillage.fit(["Reduced Till", "Till"])
    encoder_cover.fit(["Cover Crop", "No Cover Crop"])
    design.Sample_Depth = encoder_depth.transform(design.Sample_Depth)
    design.Tillage = encoder_tillage.transform(design.Tillage)
    design.Cover_Crop = encoder_cover.transform(design.Cover_Crop)
    design = design.drop("Fertility_Source", axis=1)
    design = design.join(fertility_label)
    chosed = design[(design.Sample_Depth == depth)].index.tolist()
    design = design.iloc[chosed, :]
    if random == True:
        design['random'] = np.random.rand(len(chosed)).tolist()

    # design = design.drop(chosed, axis = 0)
    if domain == 'prokaryotes':
        ra = pd.read_csv("/home/moyujie/lab102/Microbiome_V2/data/ASVs_prokaryotes.csv", index_col=0)
        taxa = pd.read_csv("/home/moyujie/lab102/Microbiome_V2/data/taxa_prokaryotes.csv", index_col=0)
        ra = ra.T
        ra = ra.drop('FST_78P', axis=0)

    else:
        ra = pd.read_csv("/home/moyujie/lab102/Microbiome_V2/data/ASVs_fungi.csv", index_col=0)
        taxa = pd.read_csv("/home/moyujie/lab102/Microbiome_V2/data/taxa_fungi.csv", index_col=0)
        ra = ra.T
        ra = ra.drop('FST_78F', axis=0)
    ra = ra.reset_index(drop=True)
    ra = ra.div(ra.sum(axis=1), axis=0)
    ra = ra.iloc[chosed, :]
    num = (ra == 0).astype(int).sum(axis=0)
    a = np.array(np.where(np.array(num) <= (24 - time))).reshape(-1)
    ra = ra.iloc[:, a]
    taxa = taxa.iloc[a, :]

    return design,ra, taxa




def explained_variance_z_score(d_origi, d_reduced):
    z_score_d_origi = (d_origi - d_origi.mean()) / d_origi.std()

    z_score_d_reduced = (d_reduced - d_reduced.mean()) / d_reduced.std()
    score = explained_variance_score(z_score_d_origi, z_score_d_reduced)

    return score
