"""
Calc.Group_Summary.py
==========================
Calculate the following domain average values for each group:
    - Precipitation
    - Surface Temperature
    - Sensible Heat Flux
    - Latent Heat Flux
    - Water vapor convergence between 925 hPa to 850 hPa
"""

import numpy as np
import xarray as xr
import pandas as pd
import json
import os
import sys
sys.path.append('../')
import preprocessing.Preprocessing_Get_Data as PrepGD
import preprocessing.Preprocessing as Prep

# Get configuration
Config = json.load(open(os.path.join(os.path.dirname(__file__), '../config.json')))

def Get_Data():

    Data = {}

    # ================================================================
    #
    # Read variables from extracted data
    #
    # ================================================================
    
    for Var in ['PRECT', 'TS', 'LHFLX', 'SHFLX']:

        Data[Var] = {}

        for Run in ['CTL', 'DEF']:

            Data[Var][Run] = []

            for En in Config['Data_En']:

                # Print message
                print('Read data: {}, {}, En{}'.format(Var, Run, En))

                # Get data
                Data[Var][Run].append(PrepGD.Get_Simulation_Extract(
                    Run,
                    En,
                    Var,
                    Range='MC_Extended_Analysis',
                )[None, ...])
            
            # Concatenate data
            Data[Var][Run] = np.concatenate(Data[Var][Run], axis=0)

            # Calculate spatial average
            Data[Var][Run] = Prep.Calc_SpatialAverage(Data[Var][Run], LandMask='Land', RangeMask='MC_Analysis', Range_Original='MC_Extended_Analysis')

            # Calculate annual mean
            Data[Var][Run] = Prep.Calc_AnnualMean(Data[Var][Run][:, 6:-6, ...], Time_Axis=1)
    
    # ================================================================
    #
    # Read variables from MSE budget data
    #
    # ================================================================

    for Var in ['VdLvqv']:

        Data[Var] = {}

        for Run in ['CTL', 'DEF']:

            Data[Var][Run] = []

            for En in Config['Data_En']:

                # Print message
                print('Read data: {}, {}, En{}'.format(Var, Run, En))

                # Get data
                Data[Var][Run].append(xr.open_dataset(f'../output/Output_Data/MSE_Budget/MSE_Budget.1970-2005.En{En}.MC.nc')[f'{Var}_{Run}'].values[None, ...])
            
            # Concatenate data
            Data[Var][Run] = np.concatenate(Data[Var][Run], axis=0)

            # Calculate vertical integration between 925 hPa to 850 hPa
            # Note: The vertical coordinate is in order of 1000 hPa to 100 hPa. Therefore, the results have to be multiplied by -1.
            plev           = Prep.Get_Lev_Interp()
            plev_start     = np.argmin(abs(plev - 925))
            plev_end       = np.argmin(abs(plev - 850))
            Data[Var][Run] = -Prep.Calc_VI(Data[Var][Run][..., plev_start:plev_end+1], plev[plev_start:plev_end+1], -1)

    # Calculate ANO (difference between DEF and CTL) for all variables and delete unnecessary data
    for Var in Data.keys():

        Data[Var]['ANO'] = Data[Var]['DEF'] - Data[Var]['CTL']

        del Data[Var]['CTL']
        del Data[Var]['DEF']

    return Data

if (__name__ == '__main__'):

    # Get data
    Data = Get_Data()

    # Get group data
    Group = pd.read_csv('../output/Output_Data/Group/Group.1970-2005.csv')

    # Combine the group dataframe and the other variables
    for Var in Data.keys():

        Group[f'{Var}_ANO'] = Data[Var]['ANO'].ravel()

    # Output data
    Group.to_csv('../output/Output_Data/Group/Group.1970-2005.Multivariable_Summary.csv', index=False)