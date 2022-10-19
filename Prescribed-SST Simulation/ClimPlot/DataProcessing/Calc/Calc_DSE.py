# Calc_DSE.py

# Import Module
import netCDF4 as nc
import os

# Parameters
Cp = 1006
g  = 9.8

def Calc(ModelProperties):

	# Connect to path
	Data_Path = '../{ModelPath}/Data/'.format(ModelPath=ModelProperties['ModelPath'])
	Data_Path = os.path.abspath(Data_Path)

	if not os.path.exists(Data_Path + '/DSE/'):
		
		os.makedirs(Data_Path + '/DSE/')
	
	# Obtain T, Z3
	for i_Run in ['CTR', 'DEF']:
		
		for i_En in ModelProperties['Ensemble']['Default_En']:

			New_FilePath = \
			Data_Path + '/DSE/{ModelName}_{En}_DSE_{Run}.nc'.format(ModelName=ModelProperties['ModelName'], En=str(i_En).zfill(ModelProperties['Ensemble']['Default_En_Format']['ZeroFill']), Run=i_Run)

			if os.path.exists(New_FilePath):

				continue
			
			Data = {}
			for i_VarName in ['T', 'Z3']:

				RawData = nc.Dataset(Data_Path + '/{VarName}/{ModelName}_{En}_{VarName}_{Run}.nc'.format(VarName=i_VarName, ModelName=ModelProperties['ModelName'], En=str(i_En).zfill(ModelProperties['Ensemble']['Default_En_Format']['ZeroFill']), Run=i_Run))
				Data[i_VarName] = RawData.variables[i_VarName][:].squeeze()

				Var_lon  = RawData.variables['lon'][:]
				Var_lat  = RawData.variables['lat'][:]
				Var_lev  = RawData.variables['lev'][:]

				n_time   = RawData.dimensions['time'].size
				n_lon    = Var_lon.shape[0]
				n_lat    = Var_lat.shape[0]
				n_lev    = Var_lev.shape[0]
			
			# Calculate DSE
			Data['DSE'] = Cp * Data['T'] + g * Data['Z3']

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
			
			Create_NewVar               = New_File.createVariable('DSE', 'f8', ('time', 'lev', 'lat', 'lon'))
			Create_NewVar.standard_name = 'DSE'
			Create_NewVar.long_name     = 'Dry static energy'
			Create_NewVar.units         = 'J/kg'
			Create_NewVar[:]            = Data['DSE']
		
			New_File.close()

	return