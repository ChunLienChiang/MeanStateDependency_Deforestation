# Calc_PRECT.py

# Import Module
import netCDF4 as nc
import os

def Calc(ModelProperties):

	# Connect to path
	Data_Path = '../{ModelPath}/Data/'.format(ModelPath=ModelProperties['ModelPath'])
	Data_Path = os.path.abspath(Data_Path)

	if not os.path.exists(Data_Path + '/PRECT/'):
		
		os.makedirs(Data_Path + '/PRECT/')
	
	# Obtain PRECC, PRECL
	for i_Run in ['CTR', 'DEF']:
		
		for i_En in ModelProperties['Ensemble']['Default_En']:

			New_FilePath = \
			Data_Path + '/PRECT/{ModelName}_{En}_PRECT_{Run}.nc'.format(ModelName=ModelProperties['ModelName'], En=str(i_En).zfill(ModelProperties['Ensemble']['Default_En_Format']['ZeroFill']), Run=i_Run)

			if os.path.exists(New_FilePath):

				continue
			
			Data = {}
			for i_VarName in ['PRECC', 'PRECL']:

				RawData = nc.Dataset(Data_Path + '/{VarName}/{ModelName}_{En}_{VarName}_{Run}.nc'.format(VarName=i_VarName, ModelName=ModelProperties['ModelName'], En=str(i_En).zfill(ModelProperties['Ensemble']['Default_En_Format']['ZeroFill']), Run=i_Run))
				Data[i_VarName] = RawData.variables[i_VarName][:].squeeze()

				Var_lon  = RawData.variables['lon'][:]
				Var_lat  = RawData.variables['lat'][:]

				n_time   = RawData.dimensions['time'].size
				n_lon    = Var_lon.shape[0]
				n_lat    = Var_lat.shape[0]
			
			# Calculate PRECT
			Data['PRECT'] = Data['PRECC'] + Data['PRECL']

			# Create new nc file
			New_File = nc.Dataset(New_FilePath,'w', format='NETCDF4')
			
			New_File.createDimension('time', n_time)
			
			New_File.createDimension('lon', n_lon)
			New_File.createDimension('lat', n_lat)

			New_File.createVariable('lon', 'f', ('lon'))
			New_File.createVariable('lat', 'f', ('lat'))

			New_File.variables['lon'][:] = Var_lon
			New_File.variables['lat'][:] = Var_lat
			
			Create_NewVar               = New_File.createVariable('PRECT', 'f8', ('time', 'lat', 'lon'))
			Create_NewVar.standard_name = 'PRECT'
			Create_NewVar.long_name     = 'Precipitation (Total)'
			Create_NewVar.units         = 'Wm-2'
			Create_NewVar[:]            = Data['PRECT']
		
			New_File.close()

	return