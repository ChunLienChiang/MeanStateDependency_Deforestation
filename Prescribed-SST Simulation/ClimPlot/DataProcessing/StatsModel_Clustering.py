# StatsModel_Clustering.py

# Import Module
import numpy as np
import pandas as pd

# Get models
def Get_ClusterModels(**kwargs):
	n_jobs     = kwargs.get('n_jobs', 3)
	n_clusters = kwargs.get('n_clusters', 2)

	ClusterModels_Name = []
	ClusterModels_Object = []

	# Append k-Means clustering
	import sklearn.cluster
	ClusterModels_Name.append('k-Means')
	ClusterModels_Object.append(sklearn.cluster.KMeans(n_clusters=n_clusters, n_init=20, max_iter=100, tol=1e-5))
	
	# Append spectral clustering
	ClusterModels_Name.append('Spectral')
	ClusterModels_Object.append(sklearn.cluster.SpectralClustering(n_clusters=n_clusters, n_init=20, affinity='nearest_neighbors'))

	return pd.DataFrame({'Name': ClusterModels_Name, 'Object': ClusterModels_Object})

# Model calculation

def Calc_ClusterCenter(Data, Labels):
    ClusterCenter = {'Label': [], 'Center': []}

    for i_Labels in np.unique(Labels):

        Group_Mask = np.tile((Labels==i_Labels)[:, None], (1, Data.shape[-1]))
        ClusterCenter['Label'].append(i_Labels)
        ClusterCenter['Center'].append(np.nanmean(np.where(Group_Mask, Data, np.nan), axis=0))

    return ClusterCenter

def Calc_InfoEntropy(Data, Labels, **kwargs):

    TestLabels_Config = kwargs.get('TestLabels_Config', {'Normal_Min': -3, 'Normal_Max': 3, 'Normal_Interval': 1})

    # Group X by -3std to +3std of X
    # Data: ANO P

    Data_mean = np.nanmean(Data)
    Data_std  = np.nanstd(Data)

    Data_Group = np.arange(TestLabels_Config['Normal_Min'], \
                           TestLabels_Config['Normal_Max'] + TestLabels_Config['Normal_Interval'], \
                           TestLabels_Config['Normal_Interval'])

    # Calculate weighted average of Eta over each group
    InfoEntropy = 0
    for i_group in Data_Group:
        
        Group_Mask = (Data>=Data_mean+i_group*Data_std) & (Data<Data_mean+(i_group+1)*Data_std)

        # Calculate summation of Eta over each label
        InfoEta = 0
        for i_n_Labels in np.unique(Labels):
            
            Labels_Mask = (Labels==i_n_Labels)

            if (np.sum(~np.isnan(Data[Group_Mask & Labels_Mask]))==0):

                continue

            P_Labels    = np.sum(~np.isnan(Data[Group_Mask & Labels_Mask])) / np.sum(~np.isnan(Data[Group_Mask]))
            InfoEta    -= P_Labels * np.log(P_Labels)
        
        InfoEntropy += (np.sum(~np.isnan(Data[Group_Mask])) / np.sum(~np.isnan(Data))) * InfoEta

    return InfoEntropy
