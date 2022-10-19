# Calc_VP.py

# Import Module
import numpy as np
import netCDF4 as nc
import os

def Calc_VP(U, V, Lat):

	import metpy.calc
	import metpy.units
	import windspharm.tools as WPTool
	import windspharm.standard as WPStdrd
	
	uwnd, uwnd_info = WPTool.prep_data(U, 'tyx')
	vwnd, vwnd_info = WPTool.prep_data(V, 'tyx')

	lats, uwnd, vwnd = WPTool.order_latdim(Lat, uwnd, vwnd)

	vp = WPStdrd.VectorWind(uwnd, vwnd).velocitypotential()

	vp = WPTool.recover_data(vp, uwnd_info)
	
	vp = (vp - np.nanmean(vp))

	return vp

def Calc(ModelProperties):

	# Connect to path
	Data_Path = '../{ModelPath}/Data/'.format(ModelPath=ModelProperties['ModelPath'])
	Data_Path = os.path.abspath(Data_Path)

	if not os.path.exists(Data_Path + '/VP/'):
		
		os.makedirs(Data_Path + '/VP/')
	
	# Obtain U, V
	for i_Run in ['CTR', 'DEF']:
		
		for i_En in ModelProperties['Ensemble']['Default_En']:

			New_FilePath = \
			Data_Path + '/VP/{ModelName}_{En}_VP_{Run}.nc'.format(ModelName=ModelProperties['ModelName'], En=str(i_En).zfill(ModelProperties['Ensemble']['Default_En_Format']['ZeroFill']), Run=i_Run)

			if os.path.exists(New_FilePath):

				continue
			
			Data = {}

			for i_VarName in ['U', 'V']:

				RawData = nc.Dataset(Data_Path + '/{VarName}/{ModelName}_{En}_{VarName}_{Run}.nc'.format(VarName=i_VarName, ModelName=ModelProperties['ModelName'], En=str(i_En).zfill(ModelProperties['Ensemble']['Default_En_Format']['ZeroFill']), Run=i_Run))
				
				Data[i_VarName] = RawData.variables[i_VarName][:].squeeze()

				Var_lon  = RawData.variables['lon'][:]
				Var_lat  = RawData.variables['lat'][:]

				n_time   = RawData.dimensions['time'].size
				n_lon    = Var_lon.shape[0]
				n_lat    = Var_lat.shape[0]
			
			# Calculate VP
			Data['VP'] = Calc_VP(Data['U'][:, -2, :, :], Data['V'][:, -2, :, :], Var_lat)

			# Create new nc file
			New_File = nc.Dataset(New_FilePath,'w', format='NETCDF4')
			
			New_File.createDimension('time', n_time)
			
			New_File.createDimension('lon', n_lon)
			New_File.createDimension('lat', n_lat)

			New_File.createVariable('lon', 'f', ('lon'))
			New_File.createVariable('lat', 'f', ('lat'))

			New_File.variables['lon'][:] = Var_lon
			New_File.variables['lat'][:] = Var_lat
			
			Create_NewVar               = New_File.createVariable('VP200', 'f8', ('time', 'lat', 'lon'))
			Create_NewVar.standard_name = 'VP200'
			Create_NewVar.long_name     = 'Velocity potential in 200hPa'
			Create_NewVar.units         = 'm2s-1'
			Create_NewVar[:]            = Data['VP']
		
			New_File.close()

	return