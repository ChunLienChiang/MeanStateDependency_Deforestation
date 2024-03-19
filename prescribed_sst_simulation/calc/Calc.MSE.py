"""
Calc.MSE.py
==========================
Calculate moist static energy.
"""

import xarray as xr
import os
import json
import sys
sys.path.append('../')
import preprocessing.Preprocessing_Get_Data as PrepGD
import preprocessing.Preprocessing as Prep

# Get configuration
Config = json.load(open(os.path.join(os.path.dirname(__file__), '../config.json')))

# Declare constants
Cp = 1006
g  = 9.8
Lv = 2.25e6

if (__name__ == '__main__'):

	# Set output path
	File_Path = '../output/Output_Data/MSE/'
	if not (os.path.exists(File_Path)): os.makedirs(File_Path)

	for i_Run in ['CTL', 'DEF']:
	
		for i_En in Config['Data_En']:

			# Print message
			print('Read data: {}, En{}'.format(i_Run, i_En))

			# Set output file name
			File_Name_Global = 'MSE.{Run}.{Year_Start}-{Year_End}.En{En}.nc'.format(\
				Run=i_Run, \
				Year_Start=Config['Data_TimeRange'][0], \
				Year_End=Config['Data_TimeRange'][1], \
				En=i_En, \
			)
			File_Name_MC     = 'MSE.{Run}.{Year_Start}-{Year_End}.En{En}.MC_Extended.nc'.format(\
				Run=i_Run, \
				Year_Start=Config['Data_TimeRange'][0], \
				Year_End=Config['Data_TimeRange'][1], \
				En=i_En, \
			)

			if (os.path.exists(File_Path + File_Name_Global)) and (os.path.exists(File_Path + File_Name_MC)): continue

			# Get data
			T = PrepGD.Get_Simulation(i_Run, i_En, 'T')
			Z = PrepGD.Get_Simulation(i_Run, i_En, 'Z3')
			Q = PrepGD.Get_Simulation(i_Run, i_En, 'Q')

			# Calculate MSE, DSE, and Lvqv
			MSE  = Cp * T + g * Z + Lv * Q
			DSE  = Cp * T + g * Z
			Lvqv = Lv * Q

			# Output to new nc file: global
			xr.Dataset(\
				data_vars={\
					'MSE': (('time', 'lev_interp', 'lat', 'lon'), MSE), \
					'DSE': (('time', 'lev_interp', 'lat', 'lon'), DSE), \
					'Lvqv': (('time', 'lev_interp', 'lat', 'lon'), Lvqv), \
				}, \
				coords={\
					'time': PrepGD.Get_Simulation(i_Run, None, 'time'), \
					'lev_interp': Prep.Get_Lev_Interp(), \
					'lat': Prep.Get_RefData(Var='lat'), \
					'lon': Prep.Get_RefData(Var='lon'), \
				}
			).to_netcdf(File_Path + File_Name_Global)

			# Output to new nc file: MC
			Lat_Crop, Lon_Crop = Prep.Crop_Lat_Lon('MC_Extended_Analysis')

			xr.Dataset(\
				data_vars={\
					'MSE': (('time', 'lev_interp', 'lat', 'lon'), Prep.Crop_Range(MSE, 'MC_Extended_Analysis')), \
					'DSE': (('time', 'lev_interp', 'lat', 'lon'), Prep.Crop_Range(DSE, 'MC_Extended_Analysis')), \
					'Lvqv': (('time', 'lev_interp', 'lat', 'lon'), Prep.Crop_Range(Lvqv, 'MC_Extended_Analysis')), \
				}, \
				coords={\
					'time': PrepGD.Get_Simulation(i_Run, None, 'time'), \
					'lev_interp': Prep.Get_Lev_Interp(), \
					'lat': Lat_Crop, \
					'lon': Lon_Crop, \
				}
			).to_netcdf(File_Path + File_Name_MC)