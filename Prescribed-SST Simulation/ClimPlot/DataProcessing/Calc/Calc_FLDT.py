# Calc_FLDT.py

# Import Module
import netCDF4 as nc
import os

def Calc(ModelProperties):

	# Connect to path
	Data_Path = '../{ModelPath}/Data/'.format(ModelPath=ModelProperties['ModelPath'])
	Data_Path = os.path.abspath(Data_Path)

	if not os.path.exists(Data_Path + '/FLDT/'):
		
		os.makedirs(Data_Path + '/FLDT/')
	
	# Obtain FLUT, FLNT
	for i_Run in ['CTR', 'DEF']:
		
		for i_En in ModelProperties['Ensemble']['Default_En']:

			New_FilePath = \
			Data_Path + '/FLDT/{ModelName}_{En}_FLDT_{Run}.nc'.format(ModelName=ModelProperties['ModelName'], En=str(i_En).zfill(ModelProperties['Ensemble']['Default_En_Format']['ZeroFill']), Run=i_Run)

			if os.path.exists(New_FilePath):

				continue
			
			Data = {}
			for i_VarName in ['FLUT', 'FLNT']:

				RawData = nc.Dataset(Data_Path + '/{VarName}/{ModelName}_{En}_{VarName}_{Run}.nc'.format(VarName=i_VarName, ModelName=ModelProperties['ModelName'], En=str(i_En).zfill(ModelProperties['Ensemble']['Default_En_Format']['ZeroFill']), Run=i_Run))
				Data[i_VarName] = RawData.variables[i_VarName][:].squeeze()

				Var_lon  = RawData.variables['lon'][:]
				Var_lat  = RawData.variables['lat'][:]

				n_time   = RawData.dimensions['time'].size
				n_lon    = Var_lon.shape[0]
				n_lat    = Var_lat.shape[0]
			
			# Calculate FLDT = FLUT - FLNT
			Data['FLDT'] = Data['FLUT'] - Data['FLNT']

			# Create new nc file
			New_File = nc.Dataset(New_FilePath,'w', format='NETCDF4')
			
			New_File.createDimension('time', n_time)
			
			New_File.createDimension('lon', n_lon)
			New_File.createDimension('lat', n_lat)

			New_File.createVariable('lon', 'f', ('lon'))
			New_File.createVariable('lat', 'f', ('lat'))

			New_File.variables['lon'][:] = Var_lon
			New_File.variables['lat'][:] = Var_lat
			
			Create_NewVar               = New_File.createVariable('FLDT', 'f8', ('time', 'lat', 'lon'))
			Create_NewVar.standard_name = 'FLDT'
			Create_NewVar.long_name     = 'Downwelling longwave flux at top of model'
			Create_NewVar.units         = 'Wm-2'
			Create_NewVar[:]            = Data['FLDT']
		
			New_File.close()

	return