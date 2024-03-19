"""
tool_ModelOutput_Extract.py
==========================
Extract specific variables and time range from model output nc files.
"""

import numpy as np
import xarray as xr
import pandas as pd
import json
import os
import sys
sys.path.append('../../')
import preprocessing.Preprocessing_Get_Data as PrepGD
import preprocessing.Preprocessing as Prep

# Get configuration
Config = json.load(open(os.path.join(os.path.dirname(__file__), '../../config.json')))

if (__name__ == '__main__'):

	# Set output path
	File_Path = '../../output/Output_Data/ModelOutput_Extract/'
	if not (os.path.exists(File_Path)): os.makedirs(File_Path)

	# Read variable list
	VarList = pd.read_csv('VarList.csv', sep=', ', engine='python')

	for _, i_VarList in VarList.iterrows():

		for i_Run in ['CTL', 'DEF']:

			for i_En in Config['Data_En']:

				# Print message
				print('Read data: {}, {}, En{}'.format(i_VarList['Var'], i_Run, i_En))

				# Set output file name
				File_Name = 'ModelOutput_Extract.{Component}.{Var}.{Run}.{Year_Start}-{Year_End}.En{En}.{Range}.nc'.format(\
					Component=i_VarList['Component'], \
					Var=i_VarList['Var'], \
					Run=i_Run, \
					Year_Start=Config['Data_TimeRange'][0], \
					Year_End=Config['Data_TimeRange'][1], \
					En=i_En, \
					Range=i_VarList['Range'], \
				)

				if (os.path.exists(File_Path + File_Name)): continue

				# Get data
				Data = PrepGD.Get_Simulation(i_Run, i_En, i_VarList['Var'])

				# Crop range
				if (i_VarList['Range'] != 'Global_Analysis'):

					Data = Prep.Crop_Range(Data, i_VarList['Range'])
					Lat, Lon = Prep.Crop_Lat_Lon(i_VarList['Range'])
				
				# Output to new nc file
				if (np.ndim(Data) == 4):

					xr.Dataset(\
						data_vars={\
							i_VarList['Var']: (('time', 'lev_interp', 'lat', 'lon'), Data), \
						}, \
						coords={\
							'time': PrepGD.Get_Simulation(i_Run, None, 'time'), \
							'lev_interp': Prep.Get_Lev_Interp(), \
							'lat': Lat, \
							'lon': Lon, \
						}
					).to_netcdf(File_Path + File_Name)
				
				elif (np.ndim(Data) == 3):
					
					xr.Dataset(\
						data_vars={\
							i_VarList['Var']: (('time', 'lat', 'lon'), Data), \
						}, \
						coords={\
							'time': PrepGD.Get_Simulation(i_Run, None, 'time'), \
							'lat': Lat, \
							'lon': Lon, \
						}
					).to_netcdf(File_Path + File_Name)