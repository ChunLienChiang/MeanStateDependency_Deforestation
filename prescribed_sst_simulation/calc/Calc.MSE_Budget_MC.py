"""
Calc.MSE_Budget_MC.py
==========================
Calculate moist static energy budget in the MC region.
"""

import numpy as np
import scipy as sp
import xarray as xr
import metpy.calc as mpcalc
import os
import json
import warnings
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

	warnings.simplefilter("ignore", category=RuntimeWarning)

	# Set output path
	Data_File_Path = '../output/Output_Data/MSE/'
	File_Path      = '../output/Output_Data/MSE_Budget/'
	if not (os.path.exists(File_Path)): os.makedirs(File_Path)

	# Get latitude and longitude, and calculate dx and dy
	Lat_Crop, Lon_Crop = Prep.Crop_Lat_Lon('MC_Extended_Analysis')
	dx, dy = mpcalc.lat_lon_grid_deltas(Lon_Crop, Lat_Crop)
	dx, dy = dx[None, None, ...], dy[None, None, ...]
	dz     = np.diff(Prep.Get_Lev_Interp(), n=1) * 100

	for i_En in Config['Data_En']:

		# Print message
		print('Read data: En{}'.format(i_En))

		File_Name = 'MSE_Budget.{Year_Start}-{Year_End}.En{En}.MC.nc'.format(\
			Year_Start=Config['Data_TimeRange'][0], \
			Year_End=Config['Data_TimeRange'][1], \
			En=i_En, \
		)

		if (os.path.exists(File_Path + File_Name)): continue

		"""
		=========================================================================
		Get Data
		From output data and original dataset.
		=========================================================================
		"""

		U       = {}
		V       = {}
		OMEGA   = {}
		MSE     = {}
		DSE     = {}
		Lvqv    = {}
		
		for i_Run in ['CTL', 'DEF']:

			# Set output file name
			Data_File_Name = 'MSE.{Run}.{Year_Start}-{Year_End}.En{En}.MC_Extended.nc'.format(\
				Run=i_Run, \
				Year_Start=Config['Data_TimeRange'][0], \
				Year_End=Config['Data_TimeRange'][1], \
				En=i_En, \
			)
			
			if not (os.path.exists(Data_File_Path + Data_File_Name)): raise ValueError('The MSE output data file is not existed: {}'.format(Data_File_Path + Data_File_Name))
	
			# Get data: U, V, and Omega
			U[i_Run]     = PrepGD.Get_Simulation_Extract(i_Run, i_En, 'U', 'MC_Extended_Analysis')
			V[i_Run]     = PrepGD.Get_Simulation_Extract(i_Run, i_En, 'V', 'MC_Extended_Analysis')
			OMEGA[i_Run] = PrepGD.Get_Simulation_Extract(i_Run, i_En, 'OMEGA', 'MC_Extended_Analysis')
			
			U[i_Run]     = Prep.Calc_AnnualMean(U[i_Run][6:-6, ...], 0)
			V[i_Run]     = Prep.Calc_AnnualMean(V[i_Run][6:-6, ...], 0)
			OMEGA[i_Run] = Prep.Calc_AnnualMean(OMEGA[i_Run][6:-6, ...], 0)

			# Get data: MSE, DSE, and Lvqv
			RawData      = xr.open_dataset(Data_File_Path + Data_File_Name)
			MSE[i_Run]   = RawData['MSE'].to_numpy()
			DSE[i_Run]   = RawData['DSE'].to_numpy()
			Lvqv[i_Run]  = RawData['Lvqv'].to_numpy()

			MSE[i_Run]   = Prep.Calc_AnnualMean(MSE[i_Run][6:-6, ...], 0)
			DSE[i_Run]   = Prep.Calc_AnnualMean(DSE[i_Run][6:-6, ...], 0)
			Lvqv[i_Run]  = Prep.Calc_AnnualMean(Lvqv[i_Run][6:-6, ...], 0)
			
			Time         = np.arange(*Config['Data_TimeRange'])

		U['ANO']     = U['DEF'] - U['CTL']
		V['ANO']     = V['DEF'] - V['CTL']
		OMEGA['ANO'] = OMEGA['DEF'] - OMEGA['CTL']
		MSE['ANO']   = MSE['DEF'] - MSE['CTL']
		DSE['ANO']   = DSE['DEF'] - DSE['CTL']
		Lvqv['ANO']  = Lvqv['DEF'] - Lvqv['CTL']

		"""
		=========================================================================
		Calculate budget
		By metpy calculate the following terms of MSE, DSE, and Lvqv:
			1. negative vertical advection
			2. negative horizontal advection
			3. negative vertical flux convergence
			4. negative horizontal flux convergence
		=========================================================================
		"""

		Data = {}

		# Include MSE, DSE, and Lvqv
		Data['MSE']                         = {key: MSE[key] for key in ['CTL', 'DEF', 'ANO']}
		Data['DSE']                         = {key: DSE[key] for key in ['CTL', 'DEF', 'ANO']}
		Data['Lvqv']                        = {key: Lvqv[key] for key in ['CTL', 'DEF', 'ANO']}
		Data['Omega']                       = {key: OMEGA[key] for key in ['CTL', 'DEF', 'ANO']}

		# Calculate dV
		Data['dV']                          = {}
		Data['dV']['CTL']                   = mpcalc.divergence(U['CTL'], V['CTL'], dx=dx, dy=dy)
		Data['dV']['DEF']                   = mpcalc.divergence(U['DEF'], V['DEF'], dx=dx, dy=dy)
		Data['dV']['ANO']                   = Data['dV']['DEF'] - Data['dV']['CTL']

		# Calculate dOmega
		Data['dOmega']                      = {}
		Data['dOmega']['CTL']               = mpcalc.first_derivative(OMEGA['CTL'], axis=-3, delta=dz)
		Data['dOmega']['DEF']               = mpcalc.first_derivative(OMEGA['DEF'], axis=-3, delta=dz)
		Data['dOmega']['ANO']               = Data['dOmega']['DEF'] - Data['dOmega']['CTL']
		
		# Calculate Vdsd
		Data['Vdsd']                        = {}
		Data['Vdsd']['CTL']                 = -mpcalc.advection(DSE['CTL'], u=U['CTL'], v=V['CTL'], dx=dx, dy=dy, vertical_dim=None)
		Data['Vdsd']['DEF']                 = -mpcalc.advection(DSE['DEF'], u=U['DEF'], v=V['DEF'], dx=dx, dy=dy, vertical_dim=None)
		Data['Vdsd']['ANO']                 = Data['Vdsd']['DEF'] - Data['Vdsd']['CTL']
		Data['Vdsd']['ANO_ThDyn']           = -mpcalc.advection(DSE['ANO'], u=U['CTL'], v=V['CTL'], dx=dx, dy=dy, vertical_dim=None)
		Data['Vdsd']['ANO_Dyn']             = -mpcalc.advection(DSE['CTL'], u=U['ANO'], v=V['ANO'], dx=dx, dy=dy, vertical_dim=None)
		Data['Vdsd']['ANO_NL']              = -mpcalc.advection(DSE['ANO'], u=U['ANO'], v=V['ANO'], dx=dx, dy=dy, vertical_dim=None)
	
		# Calculate Omegadsd
		Data['Omegadsd']                    = {}
		Data['Omegadsd']['CTL']             = -mpcalc.advection(DSE['CTL'], w=OMEGA['CTL'], dz=dz)
		Data['Omegadsd']['DEF']             = -mpcalc.advection(DSE['DEF'], w=OMEGA['DEF'], dz=dz)
		Data['Omegadsd']['ANO']             = Data['Omegadsd']['DEF'] - Data['Omegadsd']['CTL']
		Data['Omegadsd']['ANO_ThDyn']       = -mpcalc.advection(DSE['ANO'], w=OMEGA['CTL'], dz=dz)
		Data['Omegadsd']['ANO_Dyn']         = -mpcalc.advection(DSE['CTL'], w=OMEGA['ANO'], dz=dz)
		Data['Omegadsd']['ANO_NL']          = -mpcalc.advection(DSE['ANO'], w=OMEGA['ANO'], dz=dz)
		
		# Calculate VdLvqv
		Data['VdLvqv']                      = {}
		Data['VdLvqv']['CTL']               = -mpcalc.advection(Lvqv['CTL'], u=U['CTL'], v=V['CTL'], dx=dx, dy=dy)
		Data['VdLvqv']['DEF']               = -mpcalc.advection(Lvqv['DEF'], u=U['DEF'], v=V['DEF'], dx=dx, dy=dy)
		Data['VdLvqv']['ANO']               = Data['VdLvqv']['DEF'] - Data['VdLvqv']['CTL']
		Data['VdLvqv']['ANO_ThDyn']         = -mpcalc.advection(Lvqv['ANO'], u=U['CTL'], v=V['CTL'], dx=dx, dy=dy)
		Data['VdLvqv']['ANO_Dyn']           = -mpcalc.advection(Lvqv['CTL'], u=U['ANO'], v=V['ANO'], dx=dx, dy=dy)
		Data['VdLvqv']['ANO_NL']            = -mpcalc.advection(Lvqv['ANO'], u=U['ANO'], v=V['ANO'], dx=dx, dy=dy)

		# Calculate OmegadLvqv
		Data['OmegadLvqv']                  = {}
		Data['OmegadLvqv']['CTL']           = -mpcalc.advection(Lvqv['CTL'], w=OMEGA['CTL'], dz=dz)
		Data['OmegadLvqv']['DEF']           = -mpcalc.advection(Lvqv['DEF'], w=OMEGA['DEF'], dz=dz)
		Data['OmegadLvqv']['ANO']           = Data['OmegadLvqv']['DEF'] - Data['OmegadLvqv']['CTL']
		Data['OmegadLvqv']['ANO_ThDyn']     = -mpcalc.advection(Lvqv['ANO'], w=OMEGA['CTL'], dz=dz)
		Data['OmegadLvqv']['ANO_Dyn']       = -mpcalc.advection(Lvqv['CTL'], w=OMEGA['ANO'], dz=dz)
		Data['OmegadLvqv']['ANO_NL']        = -mpcalc.advection(Lvqv['ANO'], w=OMEGA['ANO'], dz=dz)

		# Calculate Vdhm
		Data['Vdhm']                        = {}
		Data['Vdhm']['CTL']                 = -mpcalc.advection(MSE['CTL'], u=U['CTL'], v=V['CTL'], dx=dx, dy=dy)
		Data['Vdhm']['DEF']                 = -mpcalc.advection(MSE['DEF'], u=U['DEF'], v=V['DEF'], dx=dx, dy=dy)
		Data['Vdhm']['ANO']                 = Data['Vdhm']['DEF'] - Data['Vdhm']['CTL']
		Data['Vdhm']['ANO_ThDyn']           = -mpcalc.advection(MSE['ANO'], u=U['CTL'], v=V['CTL'], dx=dx, dy=dy)
		Data['Vdhm']['ANO_Dyn']             = -mpcalc.advection(MSE['CTL'], u=U['ANO'], v=V['ANO'], dx=dx, dy=dy)
		Data['Vdhm']['ANO_NL']              = -mpcalc.advection(MSE['ANO'], u=U['ANO'], v=V['ANO'], dx=dx, dy=dy)
		
		# Calculate Omegadhm
		Data['Omegadhm']                    = {}
		Data['Omegadhm']['CTL']             = -mpcalc.advection(MSE['CTL'], w=OMEGA['CTL'], dz=dz)
		Data['Omegadhm']['DEF']             = -mpcalc.advection(MSE['DEF'], w=OMEGA['DEF'], dz=dz)
		Data['Omegadhm']['ANO']             = Data['Omegadhm']['DEF'] - Data['Omegadhm']['CTL']
		Data['Omegadhm']['ANO_ThDyn']       = -mpcalc.advection(MSE['ANO'], w=OMEGA['CTL'], dz=dz)
		Data['Omegadhm']['ANO_Dyn']         = -mpcalc.advection(MSE['CTL'], w=OMEGA['ANO'], dz=dz)
		Data['Omegadhm']['ANO_NL']          = -mpcalc.advection(MSE['ANO'], w=OMEGA['ANO'], dz=dz)
		
		# ==========================================================================================================
		#
		# Output to new nc file
		#
		# ==========================================================================================================

		Output_2d_Data = {}
		Output_VI_Data = {}
		
		for i_Var in list(Data.keys()):

			for i_Run in ['CTL', 'DEF', 'ANO', 'ANO_ThDyn', 'ANO_Dyn', 'ANO_NL']:
				
				if not (i_Run in Data[i_Var]):

					continue

				Data_Spatial_Average = Prep.Calc_SpatialAverage(Data[i_Var][i_Run], LandMask='Land', RangeMask='MC_Analysis', Range_Original='MC_Extended_Analysis')
				
				Output_2d_Data['{}_{}'.format(i_Var, i_Run)] = xr.DataArray(\
					dims=('time', 'lev_interp'), \
					data=Data_Spatial_Average, \
				)

				Output_VI_Data['{}_{}'.format(i_Var, i_Run)] = xr.DataArray(\
					dims=('time'),
					data=sp.integrate.simpson(
						Data_Spatial_Average[:, :-3][:, ::-1],
						x=Prep.Get_Lev_Interp()[:-3][::-1] * 100,
						axis=-1,
					) * (1/g),
				)
		
		xr.Dataset(
			data_vars=Output_2d_Data, 
			coords={
				'time': Time, 
				'lev_interp': Prep.Get_Lev_Interp(),
			}
		).to_netcdf(File_Path + File_Name)

		xr.Dataset(
			data_vars=Output_VI_Data, 
			coords={
				'time': Time, 
			},
			attrs={
				'description': 'The vertical integral of the spatial average of the variable. The vertical integral is only calculated from the surface to the 70hPa (~375K potential temperature).',
			}
		).to_netcdf(File_Path + File_Name.replace('.nc', '.VI.nc'))