# File_MSEBudget.py

# Import Module
import numpy as np
import pandas as pd

import sys

sys.path.append('..')
import ClimPlot.DataProcessing.DataProcessing as DataProc
import ClimPlot.GetData.GetData as GetData

def Get_Data_Config(ModelProperties, **kwargs):

	Data_Config = {}
	Data_Config['Year_Start'] = 1970
	Data_Config['Year_End']   = 2004

	Data_Config['n_time']       = (Data_Config['Year_End'] - Data_Config['Year_Start'] + 1) * 12
	Data_Config['n_year']       = (Data_Config['Year_End'] - Data_Config['Year_Start'] + 1)

	Data_Config['Lon'], \
	Data_Config['Lat']          = DataProc.Get_LonLat(ModelProperties, Range=Main_Data_Range)
	Data_Config['Lev_Interp']   = DataProc.Get_Lev_Interp()

	Data_Config['n_Lon'], \
	Data_Config['n_Lat']        = Data_Config['Lon'].shape[0], Data_Config['Lat'].shape[0]
	Data_Config['n_Lev_Interp'] = len(Data_Config['Lev_Interp'])

	if ('Ensemble' in ModelProperties):

		Data_Config['n_En']       = len(ModelProperties['Ensemble']['Default_En'])
		Data_Config['Dim_Time']   = 1

	else:

		Data_Config['Dim_Time']   = 0

	return Data_Config

def Get_ModelProperties(Model):

	import json
	Json_File            = open('ModelProperties.json')
	ModelProperties      = json.load(Json_File)['Model'][Model]

	return ModelProperties

def Get_MergeData(ModelProperties):

	Table_GetMergeData = pd.DataFrame({'VarName': [], 'VarName_File': [], 'Range': []})

	for i_Var in ['MSE', 'DSE', 'T', 'Q', 'U', 'V', 'OMEGA', 'FLNS', 'FLNT', 'FSNS', 'FSNT', 'LHFLX', 'SHFLX']:
		
		Table_GetMergeData = pd.concat([Table_GetMergeData, pd.DataFrame({'VarName': [i_Var], 'VarName_File': [i_Var], 'Range': [Main_Data_Range]})])

	Data = {}

	for _, i_Var in Table_GetMergeData.iterrows():

		print(i_Var['VarName'])

		Var_Config = {'VarName': i_Var['VarName'], 'VarName_File': i_Var['VarName_File'], 'Range': i_Var['Range']}

		# Set time range
		Var_Config['Data_TimeRange'] = [Data_Config['Year_Start'], Data_Config['Year_End']]
		
		if ('YearShift' in ModelProperties['Properties']):
			
			Var_Config['YearShift'] = ModelProperties['Properties']['YearShift']

		Data[i_Var['VarName']] = GetData.GetMergeData(ModelProperties, Var_Config)
		Data[i_Var['VarName']]['ANO'] = Data[i_Var['VarName']]['DEF'] - Data[i_Var['VarName']]['CTR']

	return Data

def Calc_MSEBudget(Data):

	Cp = 1004
	Lv = 2.26e6

	import xarray as xr
	import metpy.calc as mpcalc

	dx, dy = mpcalc.lat_lon_grid_deltas(Data_Config['Lon'], Data_Config['Lat'])
	dx, dy = dx[None, None, None, ...], dy[None, None, None, ...]
	dz     = np.diff(Data_Config['Lev_Interp'], n=1) * 100

	# Calculate normal Omega
	Data['NOmega']                    = {}
	Data['NOmega']['CTR']             = np.abs(Data['OMEGA']['CTR'][:, :, 1:, :, :]) / \
										(np.linalg.norm(Data['OMEGA']['CTR'][:, :, 1:, :, :]*np.sqrt(-dz)[None, None, :, None, None], ord=2, axis=-3)[:, :, None, :, :] / \
										np.sqrt((Data_Config['Lev_Interp'][1] - Data_Config['Lev_Interp'][-1])*100))
	
	Data['NOmega']['DEF']             = np.abs(Data['OMEGA']['DEF'][:, :, 1:, :, :]) / \
										(np.linalg.norm(Data['OMEGA']['DEF'][:, :, 1:, :, :]*np.sqrt(-dz)[None, None, :, None, None], ord=2, axis=-3)[:, :, None, :, :] / \
										np.sqrt((Data_Config['Lev_Interp'][1] - Data_Config['Lev_Interp'][-1])*100))
	Data['NOmega']['ANO']             = Data['NOmega']['DEF'] - Data['NOmega']['CTR']
	Dim_Tuple     = list(Data['NOmega']['CTR'].shape)
	Dim_Tuple[-3] = 1
	Data['NOmega']['CTR']             = np.concatenate((np.full(Dim_Tuple, np.nan), Data['NOmega']['CTR']), axis=-3)
	Data['NOmega']['DEF']             = np.concatenate((np.full(Dim_Tuple, np.nan), Data['NOmega']['DEF']), axis=-3)
	Data['NOmega']['ANO']             = np.concatenate((np.full(Dim_Tuple, np.nan), Data['NOmega']['ANO']), axis=-3)
	
	Data_xarray = {}

	for i_VarName in ['MSE', 'DSE', 'T', 'Q', 'U', 'V', 'OMEGA', 'NOmega']:

		Data_xarray[i_VarName] = {}

		for i_Run in Data[i_VarName]:

			Data_xarray[i_VarName][i_Run] = xr.DataArray(data=Data[i_VarName][i_Run], dims=['En', 'Time', 'Lev', 'Lat', 'Lon'])
	
	# Calculate dV
	Data['dV']                        = {}
	Data['dV']['CTR']                 = mpcalc.divergence(Data_xarray['U']['CTR'], Data_xarray['V']['CTR'], dx=dx, dy=dy).values
	Data['dV']['DEF']                 = mpcalc.divergence(Data_xarray['U']['DEF'], Data_xarray['V']['DEF'], dx=dx, dy=dy).values
	Data['dV']['ANO']                 = Data['dV']['DEF'] - Data['dV']['CTR']

	# Calculate dOmega
	Data['dOmega']                    = {}
	Data['dOmega']['CTR']             = mpcalc.first_derivative(Data_xarray['OMEGA']['CTR'], axis=-3, delta=dz).magnitude
	Data['dOmega']['DEF']             = mpcalc.first_derivative(Data_xarray['OMEGA']['DEF'], axis=-3, delta=dz).magnitude
	Data['dOmega']['ANO']             = Data['dOmega']['DEF'] - Data['dOmega']['CTR']
	
	# Calculate VdT
	Data['VdT']                       = {}
	Data['VdT']['CTR']                = -mpcalc.advection(Data_xarray['T']['CTR'], u=Data_xarray['U']['CTR'], v=Data_xarray['V']['CTR'], dx=dx, dy=dy).values * Cp
	Data['VdT']['DEF']                = -mpcalc.advection(Data_xarray['T']['DEF'], u=Data_xarray['U']['DEF'], v=Data_xarray['V']['DEF'], dx=dx, dy=dy).values * Cp
	Data['VdT']['ANO']                = Data['VdT']['DEF'] - Data['VdT']['CTR']
	Data['VdT']['ANO_T']              = -mpcalc.advection(Data_xarray['T']['ANO'], u=Data_xarray['U']['CTR'], v=Data_xarray['V']['CTR'], dx=dx, dy=dy).values * Cp
	Data['VdT']['ANO_V']              = -mpcalc.advection(Data_xarray['T']['CTR'], u=Data_xarray['U']['ANO'], v=Data_xarray['V']['ANO'], dx=dx, dy=dy).values * Cp
	Data['VdT']['ANO_NonLinear']      = -mpcalc.advection(Data_xarray['T']['ANO'], u=Data_xarray['U']['ANO'], v=Data_xarray['V']['ANO'], dx=dx, dy=dy).values * Cp

	# Calculate OmegadT
	Data['OmegadT']                   = {}
	Data['OmegadT']['CTR']            = -mpcalc.advection(Data_xarray['T']['CTR'], w=Data_xarray['OMEGA']['CTR'], dz=dz).values * Cp
	Data['OmegadT']['DEF']            = -mpcalc.advection(Data_xarray['T']['DEF'], w=Data_xarray['OMEGA']['DEF'], dz=dz).values * Cp
	Data['OmegadT']['ANO']            = Data['OmegadT']['DEF'] - Data['OmegadT']['CTR']
	Data['OmegadT']['ANO_T']          = -mpcalc.advection(Data_xarray['T']['ANO'], w=Data_xarray['OMEGA']['CTR'], dz=dz).values * Cp
	Data['OmegadT']['ANO_Omega']      = -mpcalc.advection(Data_xarray['T']['CTR'], w=Data_xarray['OMEGA']['ANO'], dz=dz).values * Cp
	Data['OmegadT']['ANO_NonLinear']  = -mpcalc.advection(Data_xarray['T']['ANO'], w=Data_xarray['OMEGA']['ANO'], dz=dz).values * Cp

	# Calculate VdQ
	Data['VdQ']                       = {}
	Data['VdQ']['CTR']                = -mpcalc.advection(Data_xarray['Q']['CTR'], u=Data_xarray['U']['CTR'], v=Data_xarray['V']['CTR'], dx=dx, dy=dy).values * Lv
	Data['VdQ']['DEF']                = -mpcalc.advection(Data_xarray['Q']['DEF'], u=Data_xarray['U']['DEF'], v=Data_xarray['V']['DEF'], dx=dx, dy=dy).values * Lv
	Data['VdQ']['ANO']                = Data['VdQ']['DEF'] - Data['VdQ']['CTR']
	Data['VdQ']['ANO_Q']              = -mpcalc.advection(Data_xarray['Q']['ANO'], u=Data_xarray['U']['CTR'], v=Data_xarray['V']['CTR'], dx=dx, dy=dy).values * Lv
	Data['VdQ']['ANO_V']              = -mpcalc.advection(Data_xarray['Q']['CTR'], u=Data_xarray['U']['ANO'], v=Data_xarray['V']['ANO'], dx=dx, dy=dy).values * Lv
	Data['VdQ']['ANO_NonLinear']      = -mpcalc.advection(Data_xarray['Q']['ANO'], u=Data_xarray['U']['ANO'], v=Data_xarray['V']['ANO'], dx=dx, dy=dy).values * Lv

	# Calculate OmegadQ
	Data['OmegadQ']                   = {}
	Data['OmegadQ']['CTR']            = -mpcalc.advection(Data_xarray['Q']['CTR'], w=Data_xarray['OMEGA']['CTR'], dz=dz).values * Lv
	Data['OmegadQ']['DEF']            = -mpcalc.advection(Data_xarray['Q']['DEF'], w=Data_xarray['OMEGA']['DEF'], dz=dz).values * Lv
	Data['OmegadQ']['ANO']            = Data['OmegadQ']['DEF'] - Data['OmegadQ']['CTR']
	Data['OmegadQ']['ANO_Q']          = -mpcalc.advection(Data_xarray['Q']['ANO'], w=Data_xarray['OMEGA']['CTR'], dz=dz).values * Lv
	Data['OmegadQ']['ANO_Omega']      = -mpcalc.advection(Data_xarray['Q']['CTR'], w=Data_xarray['OMEGA']['ANO'], dz=dz).values * Lv
	Data['OmegadQ']['ANO_NonLinear']  = -mpcalc.advection(Data_xarray['Q']['ANO'], w=Data_xarray['OMEGA']['ANO'], dz=dz).values * Lv

	# Calculate Vdsd
	Data['Vdsd']                      = {}
	Data['Vdsd']['CTR']               = -mpcalc.advection(Data_xarray['DSE']['CTR'], u=Data_xarray['U']['CTR'], v=Data_xarray['V']['CTR'], dx=dx, dy=dy).values
	Data['Vdsd']['DEF']               = -mpcalc.advection(Data_xarray['DSE']['DEF'], u=Data_xarray['U']['DEF'], v=Data_xarray['V']['DEF'], dx=dx, dy=dy).values
	Data['Vdsd']['ANO']               = Data['Vdsd']['DEF'] - Data['Vdsd']['CTR']
	Data['Vdsd']['ANO_sd']            = -mpcalc.advection(Data_xarray['DSE']['ANO'], u=Data_xarray['U']['CTR'], v=Data_xarray['V']['CTR'], dx=dx, dy=dy).values
	Data['Vdsd']['ANO_V']             = -mpcalc.advection(Data_xarray['DSE']['CTR'], u=Data_xarray['U']['ANO'], v=Data_xarray['V']['ANO'], dx=dx, dy=dy).values
	Data['Vdsd']['ANO_NonLinear']     = -mpcalc.advection(Data_xarray['DSE']['ANO'], u=Data_xarray['U']['ANO'], v=Data_xarray['V']['ANO'], dx=dx, dy=dy).values
	
	# Calculate Omegadsd
	Data['Omegadsd']                  = {}
	Data['Omegadsd']['CTR']           = -mpcalc.advection(Data_xarray['DSE']['CTR'], w=Data_xarray['OMEGA']['CTR'], dz=dz).values
	Data['Omegadsd']['DEF']           = -mpcalc.advection(Data_xarray['DSE']['DEF'], w=Data_xarray['OMEGA']['DEF'], dz=dz).values
	Data['Omegadsd']['ANO']           = Data['Omegadsd']['DEF'] - Data['Omegadsd']['CTR']
	Data['Omegadsd']['ANO_sd']        = -mpcalc.advection(Data_xarray['DSE']['ANO'], w=Data_xarray['OMEGA']['CTR'], dz=dz).values
	Data['Omegadsd']['ANO_Omega']     = -mpcalc.advection(Data_xarray['DSE']['CTR'], w=Data_xarray['OMEGA']['ANO'], dz=dz).values
	Data['Omegadsd']['ANO_NonLinear'] = -mpcalc.advection(Data_xarray['DSE']['ANO'], w=Data_xarray['OMEGA']['ANO'], dz=dz).values
	
	# Calculate Vdhm
	Data['Vdhm']                      = {}
	Data['Vdhm']['CTR']               = -mpcalc.advection(Data_xarray['MSE']['CTR'], u=Data_xarray['U']['CTR'], v=Data_xarray['V']['CTR'], dx=dx, dy=dy).values
	Data['Vdhm']['DEF']               = -mpcalc.advection(Data_xarray['MSE']['DEF'], u=Data_xarray['U']['DEF'], v=Data_xarray['V']['DEF'], dx=dx, dy=dy).values
	Data['Vdhm']['ANO']               = Data['Vdhm']['DEF'] - Data['Vdhm']['CTR']
	Data['Vdhm']['ANO_hm']            = -mpcalc.advection(Data_xarray['MSE']['ANO'], u=Data_xarray['U']['CTR'], v=Data_xarray['V']['CTR'], dx=dx, dy=dy).values
	Data['Vdhm']['ANO_V']             = -mpcalc.advection(Data_xarray['MSE']['CTR'], u=Data_xarray['U']['ANO'], v=Data_xarray['V']['ANO'], dx=dx, dy=dy).values
	Data['Vdhm']['ANO_NonLinear']     = -mpcalc.advection(Data_xarray['MSE']['ANO'], u=Data_xarray['U']['ANO'], v=Data_xarray['V']['ANO'], dx=dx, dy=dy).values
	
	# Calculate Omegadhm
	Data['Omegadhm']                  = {}
	Data['Omegadhm']['CTR']           = -mpcalc.advection(Data_xarray['MSE']['CTR'], w=Data_xarray['OMEGA']['CTR'], dz=dz).values
	Data['Omegadhm']['DEF']           = -mpcalc.advection(Data_xarray['MSE']['DEF'], w=Data_xarray['OMEGA']['DEF'], dz=dz).values
	Data['Omegadhm']['ANO']           = Data['Omegadhm']['DEF'] - Data['Omegadhm']['CTR']
	Data['Omegadhm']['ANO_hm']        = -mpcalc.advection(Data_xarray['MSE']['ANO'], w=Data_xarray['OMEGA']['CTR'], dz=dz).values
	Data['Omegadhm']['ANO_Omega']     = -mpcalc.advection(Data_xarray['MSE']['CTR'], w=Data_xarray['OMEGA']['ANO'], dz=dz).values
	Data['Omegadhm']['ANO_NonLinear'] = -mpcalc.advection(Data_xarray['MSE']['ANO'], w=Data_xarray['OMEGA']['ANO'], dz=dz).values
	
	# Rename MSE to hm
	Data['hm'] = Data['MSE']
	del Data['MSE']

	# Rename MSE to hm
	Data['Sd'] = Data['DSE']
	del Data['DSE']

	# Rename OMEGA to Omega
	Data['Omega'] = Data['OMEGA']
	del Data['OMEGA']

	# Calculate NOmegadhm
	Data['NOmegadhm']                  = {}
	Data['NOmegadhm']['CTR']           = -mpcalc.advection(Data_xarray['MSE']['CTR'], w=Data_xarray['NOmega']['CTR'], dz=dz).values
	Data['NOmegadhm']['DEF']           = -mpcalc.advection(Data_xarray['MSE']['DEF'], w=Data_xarray['NOmega']['DEF'], dz=dz).values
	Data['NOmegadhm']['ANO']           = Data['NOmegadhm']['DEF'] - Data['NOmegadhm']['CTR']
	Data['NOmegadhm']['ANO_hm']        = -mpcalc.advection(Data_xarray['MSE']['ANO'], w=Data_xarray['NOmega']['CTR'], dz=dz).values
	Data['NOmegadhm']['ANO_Omega']     = -mpcalc.advection(Data_xarray['MSE']['CTR'], w=Data_xarray['NOmega']['ANO'], dz=dz).values
	Data['NOmegadhm']['ANO_NonLinear'] = -mpcalc.advection(Data_xarray['MSE']['ANO'], w=Data_xarray['NOmega']['ANO'], dz=dz).values
	
	# Calculate GMS (Lowest level will be ignored)
	Data['GMS']                        = {}
	Data['GMS']['CTR']                 = DataProc.Calc_VI(-mpcalc.advection(Data_xarray['MSE']['CTR'][:, :, 1:, :, :], w=Data_xarray['NOmega']['CTR'][:, :, 1:, :, :], dz=dz[1:]).values, \
														  np.array(Data_Config['Lev_Interp'])[1:] * 100, -3)
	Data['GMS']['DEF']                 = DataProc.Calc_VI(-mpcalc.advection(Data_xarray['MSE']['DEF'][:, :, 1:, :, :], w=Data_xarray['NOmega']['DEF'][:, :, 1:, :, :], dz=dz[1:]).values, \
														  np.array(Data_Config['Lev_Interp'])[1:] * 100, -3)
	Data['GMS']['ANO']                 = Data['GMS']['DEF'] - Data['GMS']['CTR']
	Data['GMS']['ANO_hm']              = DataProc.Calc_VI(-mpcalc.advection(Data_xarray['MSE']['ANO'][:, :, 1:, :, :], w=Data_xarray['NOmega']['CTR'][:, :, 1:, :, :], dz=dz[1:]).values, \
														  np.array(Data_Config['Lev_Interp'])[1:] * 100, -3)
	Data['GMS']['ANO_Omega']           = DataProc.Calc_VI(-mpcalc.advection(Data_xarray['MSE']['CTR'][:, :, 1:, :, :], w=Data_xarray['NOmega']['ANO'][:, :, 1:, :, :], dz=dz[1:]).values, \
														  np.array(Data_Config['Lev_Interp'])[1:] * 100, -3)
	Data['GMS']['ANO_NonLinear']       = DataProc.Calc_VI(-mpcalc.advection(Data_xarray['MSE']['ANO'][:, :, 1:, :, :], w=Data_xarray['NOmega']['ANO'][:, :, 1:, :, :], dz=dz[1:]).values, \
														  np.array(Data_Config['Lev_Interp'])[1:] * 100, -3)

	return Data

def Create_NewFile(Data):

	import netCDF4 as nc

	for i_Run in ['CTR', 'DEF', 'ANO', 'ANO_hm', 'ANO_sd', 'ANO_V', 'ANO_Omega', 'ANO_T', 'ANO_Q', 'ANO_NonLinear']:

		# Create new nc file
		New_File = nc.Dataset('Output_Result/{}/DataAnalysis_MSEBudget/MSEBudget_{}.nc'.format(ModelProperties['ModelName'], i_Run), 'w', format='NETCDF4')
		
		New_File.createDimension('En', Data_Config['n_En'])
		New_File.createDimension('Time', Data_Config['n_time'])
		New_File.createDimension('Lev', Data_Config['n_Lev_Interp'])
		New_File.createDimension('Lat', Data_Config['n_Lat'])
		New_File.createDimension('Lon', Data_Config['n_Lon'])

		for i_VarName in Data:
			
			if not (i_Run in Data[i_VarName]):

				continue
			
			# Create new variable

			if (len(Data[i_VarName][i_Run].shape)==5):

				i_VarDim = ('En', 'Time', 'Lev', 'Lat', 'Lon')

			elif (len(Data[i_VarName][i_Run].shape)==4):

				i_VarDim = ('En', 'Time', 'Lat', 'Lon')
			
			Create_NewVar               = New_File.createVariable(i_VarName, 'f8', i_VarDim)
			Create_NewVar.standard_name = i_VarName
			Create_NewVar[:]            = Data[i_VarName][i_Run]
		
		New_File.close()

	return

###########################################

if __name__ == '__main__':

	# Set data source
	Model             = 'AMIP'
	Main_Data_Range   = 'MC_data_2'
	ModelProperties   = Get_ModelProperties(Model)
	Data_Config       = Get_Data_Config(ModelProperties)

	Data              = Get_MergeData(ModelProperties)
	Data              = Calc_MSEBudget(Data)
	Create_NewFile(Data)