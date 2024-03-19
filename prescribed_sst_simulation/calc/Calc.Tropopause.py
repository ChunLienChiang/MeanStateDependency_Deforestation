import numpy as np
import xarray as xr
import pandas as pd
import json
import os
import sys
sys.path.append('../')
import preprocessing.Preprocessing as Prep

# Get configuration
config = json.load(open(os.path.join(os.path.dirname(__file__), '../config.json')))

def get_index(list_ensemble: list) -> tuple[np.ndarray, np.ndarray, np.ndarray]:

    # Read data: grouping
    df_group = pd.read_csv('../output/Output_Data/Group/Group.1970-2005.csv')

    # Select the ensemble members
    df_group = df_group[df_group['En'].isin(list(list_ensemble))]

    # Get the index
    ensemble = df_group['En'].values
    year     = df_group['Year'].values
    group    = df_group['Group'].values

    return ensemble, year, group

def read_theta_dataset(list_ensemble: list) -> tuple[np.ndarray, np.ndarray]:

    # Create empty lists to store the dataset
    list_data_theta_CTL = []
    list_data_theta_DEF = []

    for en in list_ensemble:

        # Create the template of data path and file name
        data_path = f'../output/Output_Data/MSE_Budget/'
        data_file = f'MSE_Budget.1970-2005.En{en}.MC.nc'

        # Read the dataset
        data = xr.open_dataset(data_path + data_file)

        data_theta_CTL = data.variables['DSE_CTL'].values / 1004.0
        data_theta_DEF = data.variables['DSE_DEF'].values / 1004.0

        # Extract the dataset
        list_data_theta_CTL.append(data_theta_CTL[None, ...])
        list_data_theta_DEF.append(data_theta_DEF[None, ...])
        
        # Close the dataset
        data.close()

    # Concatenate the dataset
    data_theta_CTL = np.concatenate(list_data_theta_CTL, axis=0)
    data_theta_DEF = np.concatenate(list_data_theta_DEF, axis=0)

    return data_theta_CTL, data_theta_DEF

if (__name__ == '__main__'):

    # Get the list of ensemble members
    list_ensemble = config['Data_En']

    # Get the index of data
    ensemble, year, group = get_index(list_ensemble)

    # Read the dataset
    data_theta_CTL, data_theta_DEF = read_theta_dataset(list_ensemble)
    
    # Calculate long-term mean
    data_theta_CTL = np.nanmean(data_theta_CTL, axis=(0, 1))
    data_theta_DEF = np.nanmean(data_theta_DEF, axis=(0, 1))

    print(data_theta_CTL, data_theta_DEF)