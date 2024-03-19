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

def read_mse_budget_dataset(
        list_ensemble: list,
    ) -> tuple[
        np.ndarray,
        np.ndarray,
        np.ndarray,
        np.ndarray,
    ]:

    # Create empty lists to store the dataset
    list_data_Omegadhm_CTL = []
    list_data_Vdhm_CTL     = []
    list_data_Omegadhm_ANO = []
    list_data_Vdhm_ANO     = []

    for en in list_ensemble:

        # Create the template of data path and file name
        data_path = f'../output/Output_Data/MSE_Budget/'
        data_file = f'MSE_Budget.1970-2005.En{en}.MC.VI.nc'

        # Read the dataset
        data = xr.open_dataset(data_path + data_file)

        # Extract the dataset
        data_Omegadhm_CTL = data.variables['Omegadhm_CTL'].values * (-1)
        data_Vdhm_CTL     = data.variables['Vdhm_CTL'].values * (-1)
        data_Omegadhm_ANO = data.variables['Omegadhm_ANO'].values * (-1)
        data_Vdhm_ANO     = data.variables['Vdhm_ANO'].values * (-1)

        # Extract the dataset
        list_data_Omegadhm_CTL.append(data_Omegadhm_CTL[None, ...])
        list_data_Vdhm_CTL.append(data_Vdhm_CTL[None, ...])
        list_data_Omegadhm_ANO.append(data_Omegadhm_ANO[None, ...])
        list_data_Vdhm_ANO.append(data_Vdhm_ANO[None, ...])
        
        # Close the dataset
        data.close()

    # Concatenate the dataset
    data_Omegadhm_CTL = np.concatenate(list_data_Omegadhm_CTL, axis=0).flatten()
    data_Vdhm_CTL     = np.concatenate(list_data_Vdhm_CTL, axis=0).flatten()
    data_Omegadhm_ANO = np.concatenate(list_data_Omegadhm_ANO, axis=0).flatten()
    data_Vdhm_ANO     = np.concatenate(list_data_Vdhm_ANO, axis=0).flatten()
    
    return data_Omegadhm_CTL, data_Vdhm_CTL, data_Omegadhm_ANO, data_Vdhm_ANO

def read_modeloutput_extract_dataset(list_ensemble: list) -> tuple[dict, dict, dict]:

    # Set the list of variables
    list_variable = [
        'LHFLX',
        'SHFLX',
        'FLUT',
        'FLDS',
        'FLNS',
        'FLNT',
        'FSNT',
        'FSDS',
        'FSNS',
        'SOLIN',
    ]

    dict_data_ctl = {}
    dict_data_def = {}
    dict_data_ano = {}

    for variable in list_variable:

        dict_data_ctl[variable] = []
        dict_data_def[variable] = []
        dict_data_ano[variable] = []

        for en in list_ensemble:

            # Create the template of data path and file name
            data_path = f'../output/Output_Data/ModelOutput_Extract/'
            data_file_ctl = f'ModelOutput_Extract.atm.{variable}.CTL.1970-2005.En{en}.MC_Extended_Analysis.nc'
            data_file_def = f'ModelOutput_Extract.atm.{variable}.DEF.1970-2005.En{en}.MC_Extended_Analysis.nc'

            # Read the dataset
            data_ctl = xr.open_dataset(data_path + data_file_ctl)
            data_def = xr.open_dataset(data_path + data_file_def)

            # Extract the dataset
            data_ctl = data_ctl.variables[variable].values
            data_def = data_def.variables[variable].values

            # Calculate the spatial average
            data_ctl = Prep.Calc_SpatialAverage(data_ctl, LandMask='Land', RangeMask='MC_Analysis', Range_Original='MC_Extended_Analysis')
            data_def = Prep.Calc_SpatialAverage(data_def, LandMask='Land', RangeMask='MC_Analysis', Range_Original='MC_Extended_Analysis')

            # Select enso year
            data_ctl = data_ctl[6:-6, ...]
            data_def = data_def[6:-6, ...]

            # Calculate annual mean
            data_ctl = np.nanmean(data_ctl.reshape(-1, 12, *data_ctl.shape[1:]), axis=1)
            data_def = np.nanmean(data_def.reshape(-1, 12, *data_def.shape[1:]), axis=1)

            # Calculate the anomaly
            data_ano = data_def - data_ctl

            # Append the dataset
            dict_data_ctl[variable].append(data_ctl[None, ...])
            dict_data_def[variable].append(data_def[None, ...])
            dict_data_ano[variable].append(data_ano[None, ...])

        # Concatenate the dataset
        dict_data_ctl[variable] = np.concatenate(dict_data_ctl[variable], axis=0).flatten()
        dict_data_def[variable] = np.concatenate(dict_data_def[variable], axis=0).flatten()
        dict_data_ano[variable] = np.concatenate(dict_data_ano[variable], axis=0).flatten()

    # Calculate necessary variables
    dict_data_ctl['FSUT'] = dict_data_ctl['SOLIN'] - dict_data_ctl['FSNT']
    dict_data_def['FSUT'] = dict_data_def['SOLIN'] - dict_data_def['FSNT']
    dict_data_ano['FSUT'] = dict_data_ano['SOLIN'] - dict_data_ano['FSNT']
    dict_data_ctl['FSUS'] = dict_data_ctl['FSDS'] - dict_data_ctl['FSNS']
    dict_data_def['FSUS'] = dict_data_def['FSDS'] - dict_data_def['FSNS']
    dict_data_ano['FSUS'] = dict_data_ano['FSDS'] - dict_data_ano['FSNS']
    dict_data_ctl['FLUS'] = dict_data_ctl['FLNS'] + dict_data_ctl['FLDS']
    dict_data_def['FLUS'] = dict_data_def['FLNS'] + dict_data_def['FLDS']
    dict_data_ano['FLUS'] = dict_data_ano['FLNS'] + dict_data_ano['FLDS']

    return dict_data_ctl, dict_data_def, dict_data_ano

def output_file(
        data_index: tuple[np.ndarray, np.ndarray, np.ndarray],
        output_file: str,
        data_Omegadhm: np.ndarray,
        data_Vdhm: np.ndarray,
        data_lhflx: np.ndarray,
        data_shflx: np.ndarray,
        data_netlw: np.ndarray,
        data_netsw: np.ndarray,
        data_toa_sw_downward: np.ndarray,
        data_toa_sw_upward: np.ndarray,
        data_toa_lw_upward: np.ndarray,
        data_toa_lw_downward: np.ndarray,
        data_surface_sw_downward: np.ndarray,
        data_surface_sw_upward: np.ndarray,
        data_surface_lw_upward: np.ndarray,
        data_surface_lw_downward: np.ndarray,
    ):

    # Set the output path and file name
    output_path = f'../output/Output_Data/MSE_Budget_Summary/'

    if not os.path.exists(output_path):

        os.makedirs(output_path)

    # Create dataframe
    df = pd.DataFrame({
        'En'                      : data_index[0],
        'Year'                    : data_index[1],
        'Group'                   : data_index[2],
        'vertical_mse_advection'  : data_Omegadhm,
        'horizontal_mse_advection': data_Vdhm,
        'latent_heat_flux'        : data_lhflx,
        'sensible_heat_flux'      : data_shflx
        'net_shortwave'           : data_netsw,
        'residual'                : - data_Omegadhm - data_Vdhm - data_lhflx - data_shflx - data_netlw - data_netsw,
        'toa_sw_downward'         : data_toa_sw_downward,
        'toa_sw_upward'           : data_toa_sw_upward,
        'toa_lw_upward'           : data_toa_lw_upward,
        'toa_lw_downward'         : data_toa_lw_downward,
        'surface_sw_downward'     : data_surface_sw_downward,
        'surface_sw_upward'       : data_surface_sw_upward,
        'surface_lw_upward'       : data_surface_lw_upward,
        'surface_lw_downward_'    : data_surface_lw_downward,
    })

    # Save the dataframe
    df.to_csv(output_path + output_file, index=False)

    return

if (__name__ == '__main__'):

    # Get the list of ensemble members
    list_ensemble = config['Data_En']

    # Get the index of data
    print('Get the index of data')
    ensemble, year, group = get_index(list_ensemble)

    # Read the dataset
    print('Read the dataset: MSE Budget')
    data_Omegadhm_CTL, data_Vdhm_CTL, data_Omegadhm_ANO, data_Vdhm_ANO = read_mse_budget_dataset(list_ensemble)
    
    # Read the dataset
    print('Read the dataset: Model Output Extract')
    data_modeloutput_ctl, data_modeloutput_def, data_modeloutput_ano = read_modeloutput_extract_dataset(list_ensemble)

    # Output the CTL to file
    print('Output the CTL to file')
    output_file(
        (ensemble, year, group),
        f'MSE_Budget_Summary.ctl.1970-2005.MC.csv',
        data_Omegadhm_CTL,
        data_Vdhm_CTL,
        data_modeloutput_ctl['LHFLX'],
        data_modeloutput_ctl['SHFLX'],
        data_modeloutput_ctl['FLNS'] - data_modeloutput_ctl['FLNT'],
        data_modeloutput_ctl['FSNT'] - data_modeloutput_ctl['FSNS'],
        data_modeloutput_ctl['SOLIN'],
        data_modeloutput_ctl['FSUT'],
        data_modeloutput_ctl['FLUT'],
        np.zeros_like(data_modeloutput_ctl['FLUT']),
        data_modeloutput_ctl['FSDS'],
        data_modeloutput_ctl['FSUS'],
        data_modeloutput_ctl['FLUS'],
        data_modeloutput_ctl['FLDS'],
    )

    # Output the ANO to file
    output_file(
        (ensemble, year, group),
        f'MSE_Budget_Summary.ano.1970-2005.MC.csv',
        data_Omegadhm_ANO,
        data_Vdhm_ANO,
        data_modeloutput_ano['LHFLX'],
        data_modeloutput_ano['SHFLX'],
        data_modeloutput_ano['FLNS'] - data_modeloutput_ano['FLNT'],
        data_modeloutput_ano['FSNT'] - data_modeloutput_ano['FSNS'],
        data_modeloutput_ano['SOLIN'],
        data_modeloutput_ano['FSUT'],
        data_modeloutput_ano['FLUT'],
        np.zeros_like(data_modeloutput_ano['FLUT']),
        data_modeloutput_ano['FSDS'],
        data_modeloutput_ano['FSUS'],
        data_modeloutput_ano['FLUS'],
        data_modeloutput_ano['FLDS'],
    )