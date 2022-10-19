# Calc_MSE_Budget.py
import numpy as np
import netCDF4 as nc
import metpy.calc as mpcalc
import os
import warnings
import Preprocessing as Prep
import Preprocessing_Get_Data as PrepGetData

def Calc_HAX(X, U, V, Lat, Lon):

	dx, dy = mpcalc.lat_lon_grid_deltas(Lon, Lat)

	# d(V)
	d_V                = {}
	d_V['REF']         = -mpcalc.divergence(U['REF'], V['REF'], dx=dx[None, None, ...], dy=dy[None, None, ...])
	d_V['EXP']         = -mpcalc.divergence(U['EXP'], V['EXP'], dx=dx[None, None, ...], dy=dy[None, None, ...])
	
	# d(V * X)
	d_VX               = {}
	d_VX['REF']        = -mpcalc.divergence(U['REF'] * X['REF'], V['REF'] * X['REF'], dx=dx[None, None, ...], dy=dy[None, None, ...])
	d_VX['EXP']        = -mpcalc.divergence(U['EXP'] * X['EXP'], V['EXP'] * X['EXP'], dx=dx[None, None, ...], dy=dy[None, None, ...])

	# X * d(V)
	X_d_V              = {}
	X_d_V['REF']       = -X['REF'] * mpcalc.divergence(U['REF'], V['REF'], dx=dx[None, None, ...], dy=dy[None, None, ...])
	X_d_V['EXP']       = -X['EXP'] * mpcalc.divergence(U['EXP'], V['EXP'], dx=dx[None, None, ...], dy=dy[None, None, ...])

	# V * d(X)
	V_d_X              = {}
	V_d_X['REF']       = mpcalc.advection(X['REF'], u=U['REF'], v=V['REF'], w=None, dx=dx[None, None, ...], dy=dy[None, None, ...])
	V_d_X['EXP']       = mpcalc.advection(X['EXP'], u=U['EXP'], v=V['EXP'], w=None, dx=dx[None, None, ...], dy=dy[None, None, ...])
	V_d_X['ANO_V']     = mpcalc.advection(\
							np.nanmean(X['REF'], axis=0), \
							u=np.nanmean(U['EXP'], axis=0)-np.nanmean(U['REF'], axis=0), \
							v=np.nanmean(V['EXP'], axis=0)-np.nanmean(V['REF'], axis=0), \
							w=None, dx=dx[None, ...], dy=dy[None, ...])
	V_d_X['ANO_X']     = mpcalc.advection(\
							np.nanmean(X['EXP'], axis=0)-np.nanmean(X['REF'], axis=0), \
							u=np.nanmean(U['REF'], axis=0), v=np.nanmean(V['REF'], axis=0), \
							w=None, dx=dx[None, ...], dy=dy[None, ...])
	V_d_X['ANO_NL']    = mpcalc.advection(\
							np.nanmean(X['EXP'], axis=0)-np.nanmean(X['REF'], axis=0), \
							u=np.nanmean(U['EXP'], axis=0)-np.nanmean(U['REF'], axis=0), \
							v=np.nanmean(V['EXP'], axis=0)-np.nanmean(V['REF'], axis=0), \
							w=None, dx=dx[None, ...], dy=dy[None, ...])

	return d_V, d_VX, X_d_V, V_d_X

def Calc_VAX(X, OMEGA, Lev):

	# d(OMEGA * X)
	d_OMEGAX               = {}
	d_OMEGAX['REF']        = -mpcalc.first_derivative(OMEGA['REF'] * X['REF'], axis=1, x=Lev)
	d_OMEGAX['EXP']        = -mpcalc.first_derivative(OMEGA['EXP'] * X['EXP'], axis=1, x=Lev)

	# X * d(OMEGA)
	X_d_OMEGA              = {}
	X_d_OMEGA['REF']       = -X['REF'] * mpcalc.first_derivative(OMEGA['REF'], axis=1, x=Lev)
	X_d_OMEGA['EXP']       = -X['EXP'] * mpcalc.first_derivative(OMEGA['EXP'], axis=1, x=Lev)

	# OMEGA * d(X)
	OMEGA_d_X              = {}
	OMEGA_d_X['REF']       = mpcalc.advection(X['REF'], u=None, v=None, w=OMEGA['REF'], dz=np.diff(Lev, 1)[None, :, None, None])
	OMEGA_d_X['EXP']       = mpcalc.advection(X['EXP'], u=None, v=None, w=OMEGA['EXP'], dz=np.diff(Lev, 1)[None, :, None, None])
	OMEGA_d_X['ANO_V']     = mpcalc.advection(\
								np.nanmean(X['REF'], axis=0), \
								w=np.nanmean(OMEGA['EXP'], axis=0)-np.nanmean(OMEGA['REF'], axis=0), \
								dz=np.diff(Lev, 1)[:, None, None])
	OMEGA_d_X['ANO_X']     = mpcalc.advection(\
								np.nanmean(X['EXP'], axis=0)-np.nanmean(X['REF'], axis=0), \
								w=np.nanmean(OMEGA['REF'], axis=0), \
								dz=np.diff(Lev, 1)[:, None, None])
	OMEGA_d_X['ANO_NL']    = mpcalc.advection(\
								np.nanmean(X['EXP'], axis=0)-np.nanmean(X['REF'], axis=0), \
								w=np.nanmean(OMEGA['EXP'], axis=0)-np.nanmean(OMEGA['REF'], axis=0), \
								dz=np.diff(Lev, 1)[:, None, None])
	
	return d_OMEGAX, X_d_OMEGA, OMEGA_d_X

if (__name__ == '__main__'):
	
	warnings.simplefilter('ignore', category=RuntimeWarning)
	
	Cp = 1006
	g  = 9.8
	Lv = 2.25e6
	
	Model_Properties = Prep.Get_Model_Properties()

	# Get the list of events
	Events_List = PrepGetData.Get_Events_List()
	print('Calculated events (Macro-perturbation) list: {}. {} events in total.'.format(Events_List, len(Events_List)))

	File_Path = 'Output_Data/nc_MSE_Budget/'
	if not (os.path.exists(File_Path)): os.mkdir(File_Path)

	#MiPert_List  = np.arange(5)
	MiPert_List  = np.array([0])
	n_MiPert     = np.size(MiPert_List)

	for i_Event in Events_List:

		for i_MiPert in MiPert_List:

			if ((i_MiPert is None) or (i_MiPert == 0)):
				
				Suffix_MiPert = ''

			else:
				
				Suffix_MiPert = '_e' + str(i_MiPert)

			File_Name = 'nc_MSE_Budget_piControl_{Event:04d}{Suffix_MiPert}.nc'
			File_Name_Arg = {'Event': i_Event, 'Suffix_MiPert': Suffix_MiPert}
			
			if (os.path.exists((File_Path + File_Name).format(**File_Name_Arg))): continue

			print('Calculating MSE budget: {:04d}{}'.format(i_Event, Suffix_MiPert))

			Data = {}

			for i_Var in ['T', 'Z3', 'Q', 'U', 'V', 'OMEGA']:

				Data[i_Var] = {}

				for i_Run in ['CTL', 'DEFSF']:

					# Get data
					Data[i_Var][i_Run] = PrepGetData.Get_RestartSimulation(i_Run, i_Event, i_Var, MiPert=i_MiPert)
					Data[i_Var][i_Run] = Prep.Crop_Range(Data[i_Var][i_Run], 'MC_Analysis_Extended', Range_Original='Global_Analysis')

			Data['MSE'] = {}

			Data['MSE']['CTL']   = Data['T']['CTL'] * Cp + Data['Z3']['CTL'] * g + Data['Q']['CTL'] * Lv
			#Data['MSE']['DEF']   = Data['T']['DEF'] * Cp + Data['Z3']['DEF'] * g + Data['Q']['DEF'] * Lv
			Data['MSE']['DEFSF'] = Data['T']['DEFSF'] * Cp + Data['Z3']['DEFSF'] * g + Data['Q']['DEFSF'] * Lv
			
			Data['d_V']          = {}
			Data['d_V_MSE']      = {}
			Data['MSE_d_V']      = {}
			Data['V_d_MSE']      = {}
			Data['d_OMEGA_MSE']  = {}
			Data['MSE_d_OMEGA']  = {}
			Data['OMEGA_d_MSE']  = {}

			Data['d_V_Q']        = {}
			Data['Q_d_V']        = {}
			Data['V_d_Q']        = {}
			Data['d_OMEGA_Q']    = {}
			Data['Q_d_OMEGA']    = {}
			Data['OMEGA_d_Q']    = {}
			
			for i_DEF_Suffix in ['SF']:

				# Calculate horizontal advection of MSE
				d_V, d_V_MSE, MSE_d_V, V_d_MSE = Calc_HAX(\
					{'REF': Data['MSE']['CTL'], 'EXP': Data['MSE']['DEF'+i_DEF_Suffix]}, \
					{'REF': Data['U']['CTL'], 'EXP': Data['U']['DEF'+i_DEF_Suffix]}, \
					{'REF': Data['V']['CTL'], 'EXP': Data['V']['DEF'+i_DEF_Suffix]}, \
					*Prep.Crop_Lat_Lon(Model_Properties['Lat'], Model_Properties['Lon'], 'MC_Analysis_Extended') \
				)
				
				# Calculate vertical advection of MSE
				d_OMEGA_MSE, MSE_d_OMEGA, OMEGA_d_MSE = Calc_VAX(\
					{'REF': Data['MSE']['CTL'], 'EXP': Data['MSE']['DEF'+i_DEF_Suffix]}, \
					{'REF': Data['OMEGA']['CTL'], 'EXP': Data['OMEGA']['DEF'+i_DEF_Suffix]}, \
					Model_Properties['Lev_Interp']*100 \
				)
				
				Data['d_V']['CTL']                             = d_V['REF']
				Data['d_V']['DEF'+i_DEF_Suffix]                = d_V['EXP']
				Data['d_V_MSE']['CTL']                         = d_V_MSE['REF']
				Data['d_V_MSE']['DEF'+i_DEF_Suffix]            = d_V_MSE['EXP']
				Data['MSE_d_V']['CTL']                         = MSE_d_V['REF']
				Data['MSE_d_V']['DEF'+i_DEF_Suffix]            = MSE_d_V['EXP']
				Data['V_d_MSE']['CTL']                         = V_d_MSE['REF']
				Data['V_d_MSE']['DEF'+i_DEF_Suffix]            = V_d_MSE['EXP']
				Data['V_d_MSE']['ANO']                         = V_d_MSE['EXP'] - V_d_MSE['REF']
				Data['V_d_MSE']['ANO'+i_DEF_Suffix+'_V']       = V_d_MSE['ANO_V']
				Data['V_d_MSE']['ANO'+i_DEF_Suffix+'_MSE']     = V_d_MSE['ANO_X']
				Data['V_d_MSE']['ANO'+i_DEF_Suffix+'_NL']      = V_d_MSE['ANO_NL']
					
				Data['d_OMEGA_MSE']['CTL']                     = d_OMEGA_MSE['REF']
				Data['d_OMEGA_MSE']['DEF'+i_DEF_Suffix]        = d_OMEGA_MSE['EXP']
				Data['MSE_d_OMEGA']['CTL']                     = MSE_d_OMEGA['REF']
				Data['MSE_d_OMEGA']['DEF'+i_DEF_Suffix]        = MSE_d_OMEGA['EXP']
				Data['OMEGA_d_MSE']['CTL']                     = OMEGA_d_MSE['REF']
				Data['OMEGA_d_MSE']['DEF'+i_DEF_Suffix]        = OMEGA_d_MSE['EXP']
				Data['OMEGA_d_MSE']['ANO']                     = OMEGA_d_MSE['EXP'] - OMEGA_d_MSE['REF']
				Data['OMEGA_d_MSE']['ANO'+i_DEF_Suffix+'_V']   = OMEGA_d_MSE['ANO_V']
				Data['OMEGA_d_MSE']['ANO'+i_DEF_Suffix+'_MSE'] = OMEGA_d_MSE['ANO_X']
				Data['OMEGA_d_MSE']['ANO'+i_DEF_Suffix+'_NL']  = OMEGA_d_MSE['ANO_NL']

				# Calculate horizontal advection of moisture
				d_V, d_V_Q, Q_d_V, V_d_Q = Calc_HAX(\
					{'REF': Data['Q']['CTL'], 'EXP': Data['Q']['DEF'+i_DEF_Suffix]}, \
					{'REF': Data['U']['CTL'], 'EXP': Data['U']['DEF'+i_DEF_Suffix]}, \
					{'REF': Data['V']['CTL'], 'EXP': Data['V']['DEF'+i_DEF_Suffix]}, \
					*Prep.Crop_Lat_Lon(Model_Properties['Lat'], Model_Properties['Lon'], 'MC_Analysis_Extended') \
				)
				
				# Calculate vertical advection of moisture
				d_OMEGA_Q, Q_d_OMEGA, OMEGA_d_Q = Calc_VAX(\
					{'REF': Data['Q']['CTL'], 'EXP': Data['Q']['DEF'+i_DEF_Suffix]}, \
					{'REF': Data['OMEGA']['CTL'], 'EXP': Data['OMEGA']['DEF'+i_DEF_Suffix]}, \
					Model_Properties['Lev_Interp']*100 \
				)
				
				Data['d_V_Q']['CTL']                           = d_V_Q['REF']
				Data['d_V_Q']['DEF'+i_DEF_Suffix]              = d_V_Q['EXP']
				Data['Q_d_V']['CTL']                           = Q_d_V['REF']
				Data['Q_d_V']['DEF'+i_DEF_Suffix]              = Q_d_V['EXP']
				Data['V_d_Q']['CTL']                           = V_d_Q['REF']
				Data['V_d_Q']['DEF'+i_DEF_Suffix]              = V_d_Q['EXP']
				Data['V_d_Q']['ANO']                           = V_d_Q['EXP'] - V_d_Q['REF']
				Data['V_d_Q']['ANO'+i_DEF_Suffix+'_V']         = V_d_Q['ANO_V']
				Data['V_d_Q']['ANO'+i_DEF_Suffix+'_Q']         = V_d_Q['ANO_X']
				Data['V_d_Q']['ANO'+i_DEF_Suffix+'_NL']        = V_d_Q['ANO_NL']
					
				Data['d_OMEGA_Q']['CTL']                       = d_OMEGA_Q['REF']
				Data['d_OMEGA_Q']['DEF'+i_DEF_Suffix]          = d_OMEGA_Q['EXP']
				Data['Q_d_OMEGA']['CTL']                       = Q_d_OMEGA['REF']
				Data['Q_d_OMEGA']['DEF'+i_DEF_Suffix]          = Q_d_OMEGA['EXP']
				Data['OMEGA_d_Q']['CTL']                       = OMEGA_d_Q['REF']
				Data['OMEGA_d_Q']['DEF'+i_DEF_Suffix]          = OMEGA_d_Q['EXP']
				Data['OMEGA_d_Q']['ANO']                       = OMEGA_d_Q['EXP'] - OMEGA_d_Q['REF']
				Data['OMEGA_d_Q']['ANO'+i_DEF_Suffix+'_V']     = OMEGA_d_Q['ANO_V']
				Data['OMEGA_d_Q']['ANO'+i_DEF_Suffix+'_Q']     = OMEGA_d_Q['ANO_X']
				Data['OMEGA_d_Q']['ANO'+i_DEF_Suffix+'_NL']    = OMEGA_d_Q['ANO_NL']

			# Calculate spatial average
			for i_Var in list(Data.keys()):

				for i_Run in list(Data[i_Var].keys()):
					
					Data[i_Var][i_Run] = Prep.Calc_SpatialAverage(Data[i_Var][i_Run], 'MC_Analysis_Extended', LandFraction_Type='Land', Mask_Range=False)

			# Write to nc file
			New_File = nc.Dataset((File_Path + File_Name).format(**File_Name_Arg), mode='w', format='NETCDF4')
			
			# Create dimensions
			New_File.createDimension('Time', Data['MSE']['CTL'].shape[0])
			New_File.createDimension('Lev', Model_Properties['Lev_Interp'].shape[0])
			
			for i_Run in ['CTL', 'DEFSF']:
			
				New_Var    = New_File.createVariable('CpT_{}'.format(i_Run), np.float32, ('Time', 'Lev'))
				New_Var[:] = Data['T'][i_Run] * Cp
				
				New_Var    = New_File.createVariable('gz_{}'.format(i_Run), np.float32, ('Time', 'Lev'))
				New_Var[:] = Data['Z3'][i_Run] * g
				
				New_Var    = New_File.createVariable('Lvqv_{}'.format(i_Run), np.float32, ('Time', 'Lev'))
				New_Var[:] = Data['Q'][i_Run] * Lv
				
				New_Var    = New_File.createVariable('MSE_{}'.format(i_Run), np.float32, ('Time', 'Lev'))
				New_Var[:] = Data['MSE'][i_Run]
				
				New_Var    = New_File.createVariable('OMEGA_{}'.format(i_Run), np.float32, ('Time', 'Lev'))
				New_Var[:] = Data['OMEGA'][i_Run]

			for i_Var in ['d_V', 'd_V_MSE', 'MSE_d_V', 'V_d_MSE', 'd_OMEGA_MSE', 'MSE_d_OMEGA', 'OMEGA_d_MSE', 'd_V_Q', 'Q_d_V', 'V_d_Q', 'd_OMEGA_Q', 'Q_d_OMEGA', 'OMEGA_d_Q']:

				for i_Run in list(Data[i_Var].keys()):

					New_Var    = New_File.createVariable('{}_{}'.format(i_Var, i_Run), np.float32, ('Time', 'Lev'))
					New_Var[:] = Data[i_Var][i_Run]

			New_File.close()