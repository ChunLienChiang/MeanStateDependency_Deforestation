# Calc_MSE.py

# Import Module
import netCDF4 as nc
import os

# Parameters
Cp = 1006
g  = 9.8
Lv = 2.25e6

def Calc(ModelProperties):

	# Connect to path
	Data_Path = '../{ModelPath}/Data/'.format(ModelPath=ModelProperties['ModelPath'])
	Data_Path = os.path.abspath(Data_Path)

	if not os.path.exists(Data_Path + '/MSE/'):
		
		os.makedirs(Data_Path + '/MSE/')
	
	# Obtain T, Z3, Q
	for i_Run in ['CTR', 'DEF']:
		
		for i_En in ModelProperties['Ensemble']['Default_En']:

			New_FilePath = \
			Data_Path + '/MSE/{ModelName}_{En}_MSE_{Run}.nc'.format(ModelName=ModelProperties['ModelName'], En=str(i_En).zfill(ModelProperties['Ensemble']['Default_En_Format']['ZeroFill']), Run=i_Run)

			if os.path.exists(New_FilePath):

				continue
			
			Data = {}
			for i_VarName in ['T', 'Z3', 'Q']:

				print('Read {} {}'.format(i_Run, i_VarName))

				RawData = nc.Dataset(Data_Path + '/{VarName}/{ModelName}_{En}_{VarName}_{Run}.nc'.format(VarName=i_VarName, ModelName=ModelProperties['ModelName'], En=str(i_En).zfill(ModelProperties['Ensemble']['Default_En_Format']['ZeroFill']), Run=i_Run))
				Data[i_VarName] = RawData.variables[i_VarName][:].squeeze()

				Var_lon  = RawData.variables['lon'][:]
				Var_lat  = RawData.variables['lat'][:]
				Var_lev  = RawData.variables['lev'][:]

				n_time   = RawData.dimensions['time'].size
				n_lon    = Var_lon.shape[0]
				n_lat    = Var_lat.shape[0]
				n_lev    = Var_lev.shape[0]
			
			# Calculate MSE
			print('Calculate MSE')
			Data['MSE'] = Cp * Data['T'] + g * Data['Z3'] + Lv * Data['Q']

			# Create new nc file
			New_File = nc.Dataset(New_FilePath,'w', format='NETCDF4')
			
			New_File.createDimension('time', n_time)
			
			New_File.createDimension('lon', n_lon)
			New_File.createDimension('lat', n_lat)
			New_File.createDimension('lev', n_lev)

			New_File.createVariable('lon', 'f', ('lon'))
			New_File.createVariable('lat', 'f', ('lat'))
			New_File.createVariable('lev', 'f', ('lev'))

			New_File.variables['lon'][:] = Var_lon
			New_File.variables['lat'][:] = Var_lat
			New_File.variables['lev'][:] = Var_lev
			
			Create_NewVar               = New_File.createVariable('MSE', 'f8', ('time', 'lev', 'lat', 'lon'))
			Create_NewVar.standard_name = 'MSE'
			Create_NewVar.long_name     = 'Moist static energy'
			Create_NewVar.units         = 'J/kg2'
			Create_NewVar[:]            = Data['MSE']
		
			New_File.close()

	return