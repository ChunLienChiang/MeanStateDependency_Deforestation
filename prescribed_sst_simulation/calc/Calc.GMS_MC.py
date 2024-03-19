"""
Calc.GMS_MC.py
==========================
Calculate gross moist stability in the MC region.
"""

import numpy as np
import xarray as xr
import metpy.calc as mpcalc
import os
import json
import sys
sys.path.append('../')
import preprocessing.Preprocessing_Get_Data as PrepGD
import preprocessing.Preprocessing as Prep

# Get configuration
Config = json.load(open(os.path.join(os.path.dirname(__file__), '../config.json')))

if (__name__ == '__main__'):

	# Declare constants
	g  = 9.8

	# Get interpolated level, and calculate dz (dp)
	dz     = np.diff(Prep.Get_Lev_Interp(), n=1) * 100
	P_T    = (Prep.Get_Lev_Interp()[1] - Prep.Get_Lev_Interp()[-1])*100

	for i_En in Config['Data_En']:

		# Print message
		print('Read data: En{}'.format(i_En))

		File_Path = '../output/Output_Data/GMS/'
		File_Name = 'GMS.{Year_Start}-{Year_End}.En{En}.MC_Extended.nc'.format(\
			Year_Start=Config['Data_TimeRange'][0], \
			Year_End=Config['Data_TimeRange'][1], \
			En=i_En, \
		)
		if not (os.path.exists(File_Path)): os.makedirs(File_Path)

		if (os.path.exists(File_Path + File_Name)): continue

		"""
		=========================================================================
		Get Data
		From output data and original dataset.
		=========================================================================
		"""

		MSE    = {}
		OMEGA  = {}
		NOmega = {}
		
		for i_Run in ['CTL', 'DEF']:

			# Set output file name
			Data_File_Path = '../output/Output_Data/MSE_Budget/'
			Data_File_Name = 'MSE_Budget.{Year_Start}-{Year_End}.En{En}.MC_Extended.nc'.format(\
				Year_Start=Config['Data_TimeRange'][0], \
				Year_End=Config['Data_TimeRange'][1], \
				En=i_En, \
			)
			
			if not (os.path.exists(Data_File_Path + Data_File_Name)): raise ValueError('The MSE budget data file is not existed: {}'.format(Data_File_Path + Data_File_Name))
	
			# Get data: MSE
			RawData       = xr.open_dataset(Data_File_Path + Data_File_Name)
			MSE[i_Run]    = RawData['MSE_{}'.format(i_Run)].to_numpy()

			# Set time
			Time          = np.arange(*Config['Data_TimeRange'])

			# Get data: Omega
			OMEGA[i_Run]  = PrepGD.Get_Simulation_Extract(i_Run, i_En, 'OMEGA', 'MC_Extended_Analysis')
			OMEGA[i_Run]  = Prep.Calc_SpatialAverage(OMEGA[i_Run], LandMask='Land', RangeMask='MC_Analysis', Range_Original='MC_Extended_Analysis')
			OMEGA[i_Run]  = Prep.Calc_AnnualMean(OMEGA[i_Run][6:-6, ...], 0)
			
			# Calculate normalized omega
			NOmega[i_Run] = np.abs(OMEGA[i_Run][..., 1:]) / (np.linalg.norm(OMEGA[i_Run][..., 1:]*np.sqrt(-dz)[None, :], ord=2, axis=-1) / np.sqrt(P_T))[:, None]
			NOmega[i_Run] = np.concatenate((np.full((NOmega[i_Run].shape[0], 1), np.nan), NOmega[i_Run]), axis=1)
		
		MSE['ANO']    = MSE['DEF'] - MSE['CTL']
		OMEGA['ANO']  = OMEGA['DEF'] - OMEGA['CTL']
		NOmega['ANO'] = NOmega['DEF'] - NOmega['CTL']

		"""
		=========================================================================
		Calculate GMS
		=========================================================================
		"""

		Data = {}

		# Calculate GMS
		Data['GMS']                 = {}
		Data['GMS']['CTL']          = -Prep.Calc_VI(\
			-mpcalc.advection(MSE['CTL'][..., 1:], w=NOmega['CTL'][..., 1:], u=None, v=None, dz=dz[1:], vertical_dim=-1).filled(np.nan), \
			np.array(Prep.Get_Lev_Interp())[1:] * 100, -1, \
		) * (-g/P_T)
		Data['GMS']['DEF']          = -Prep.Calc_VI(\
			-mpcalc.advection(MSE['DEF'][..., 1:], w=NOmega['DEF'][..., 1:], u=None, v=None, dz=dz[1:], vertical_dim=-1).filled(np.nan), \
			np.array(Prep.Get_Lev_Interp())[1:] * 100, -1, \
		) * (-g/P_T)
		Data['GMS']['ANO']          = Data['GMS']['DEF'] - Data['GMS']['CTL']
		Data['GMS']['ANO_ThDyn']    = -Prep.Calc_VI(\
			-mpcalc.advection(MSE['ANO'][..., 1:], w=NOmega['CTL'][..., 1:], u=None, v=None, dz=dz[1:], vertical_dim=-1).filled(np.nan), \
			np.array(Prep.Get_Lev_Interp())[1:] * 100, -1, \
		) * (-g/P_T)
		Data['GMS']['ANO_Dyn']      = -Prep.Calc_VI(\
			-mpcalc.advection(MSE['CTL'][..., 1:], w=NOmega['ANO'][..., 1:], u=None, v=None, dz=dz[1:], vertical_dim=-1).filled(np.nan), \
			np.array(Prep.Get_Lev_Interp())[1:] * 100, -1, \
		) * (-g/P_T)
		Data['GMS']['ANO_NL']       = -Prep.Calc_VI(\
			-mpcalc.advection(MSE['ANO'][..., 1:], w=NOmega['ANO'][..., 1:], u=None, v=None, dz=dz[1:], vertical_dim=-1).filled(np.nan), \
			np.array(Prep.Get_Lev_Interp())[1:] * 100, -1, \
		) * (-g/P_T)

		# Output to new nc file
		Output_Data = {}
		
		for i_Var in list(Data.keys()):

			for i_Run in ['CTL', 'DEF', 'ANO', 'ANO_ThDyn', 'ANO_Dyn', 'ANO_NL']:
				
				if not (i_Run in Data[i_Var]):

					continue
				
				Output_Data['{}_{}'.format(i_Var, i_Run)] = xr.DataArray(\
					dims=('time'), \
					data=Data[i_Var][i_Run], \
				)

		xr.Dataset(\
			data_vars=Output_Data, \
			coords={\
				'time': Time, \
			}
		).to_netcdf(File_Path + File_Name)
