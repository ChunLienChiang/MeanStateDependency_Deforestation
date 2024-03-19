import numpy as np
import pandas as pd
import xarray as xr
import scipy as sp
import json
import os
import sys
sys.path.append('../')
import preprocessing.Preprocessing as Prep

# Get configuration
Config = json.load(open(os.path.join(os.path.dirname(__file__), '../config.json')))

if (__name__ == '__main__'):

    # Read data: grouping
    Data_PRECT  = pd.read_csv('../output/Output_Data/Group/Group.1970-2005.csv')
    Mask_Group1 = (Data_PRECT['Group']==0)
    Mask_Group2 = (Data_PRECT['Group']==3)

    # Read data: MSE budget
    # Plot

    for i_Var in ['Omegadhm']:

        for i_Run in ['ANO']:

            Data = []

            try:
                
                for i_En in Config['Data_En']:

                    Data_File_Path = '../output/Output_Data/MSE_Budget/'
                    Data_File_Name = 'MSE_Budget.1970-2005.En{En}.MC.nc'.format(En=i_En)

                    if not (os.path.exists(Data_File_Path + Data_File_Name)): 
                        
                        raise ValueError('The file {} does not exist.'.format(Data_File_Name))

                    Data.append(xr.open_dataset(Data_File_Path + Data_File_Name)['{Var}_{Run}'.format(Var=i_Var, Run=i_Run)].to_numpy())

            except:

                continue

            Data = np.concatenate(tuple(Data), axis=0)

            # Change sign
            if (i_Var in ['dV', 'Vdhm', 'Omegadhm']):
                
                Data = -Data
            
            Data_VI = sp.integrate.simpson(
                Data[:, :6][:, ::-1],
                x=Prep.Get_Lev_Interp()[:6][::-1] * 100,
                axis=-1,
            )

            print(np.round(np.nanmean(Data_VI[Mask_Group1], axis=0), 3))
            print(np.round(np.nanmean(Data_VI[Mask_Group2], axis=0), 3))