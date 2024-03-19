# Preprocessing_Get_Data.py
import numpy as np
import netCDF4 as nc
import Preprocessing as Prep
import os

Cp = 1006
g  = 9.8
Lv = 2.25e6

def Calc_VerticalInterpolation(Data, PS, **kwargs):

	RawData_File = kwargs.get('RawData_File', '/work/data/CMIP6_B_piControl/TS/b.e21.B1850.f09_g17.CMIP6-piControl.001.cam.h0.TS.000101-009912.nc')

	import Ngl

	RawData_File = nc.Dataset(RawData_File)
	
	# Interpolation
	Data_Interp = Ngl.vinth2p( \
		Data.squeeze(), RawData_File.variables['hyam'][:], RawData_File.variables['hybm'][:], \
		Prep.Get_Model_Properties()['Lev_Interp'], PS.squeeze(), 1, 1000, 1, False\
	)
	
	# Mask out invalid values
	Data_Interp = np.where(Data_Interp>1e25, np.nan, Data_Interp)

	return Data_Interp

def Get_Original_TS(Year):

	Start_Year = 100
	End_Year   = 1800

	if (isinstance(Year, int)):

		RawData   = nc.Dataset('/work/data/CMIP6_B_piControl/TS/b.e21.B1850.f09_g17.CMIP6-piControl.001.cam.h0.TS.{:04d}01-{:04d}12.nc'.format(Year//100*100, Year//100*100+99), 'r')
		Data = RawData.variables['TS'][:]
		Data = Data[(Year-Year//100*100)*12:(Year-Year//100*100+1)*12, ...]

	elif (isinstance(Year, (list, np.ndarray, np.generic))):

		Year                  = np.array(Year)
		Mask_Time             = np.full((End_Year-Start_Year+100), False)
		Mask_Time[(Year-100)] = True
		Mask_Time             = np.repeat(Mask_Time, 12)

		RawData = nc.MFDataset(['/work/data/CMIP6_B_piControl/TS/b.e21.B1850.f09_g17.CMIP6-piControl.001.cam.h0.TS.{:04d}01-{:04d}12.nc'.format(i, i+99) for i in np.arange(Start_Year, End_Year+100, 100)], 'r')
		Data = RawData.variables['TS'][Mask_Time, ...]

	return Data

def Get_RestartSimulation(Run, Year, Var, **kwargs):

	MiPert                  = kwargs.get('MiPert', None)
	Compset                 = kwargs.get('Compset', 'piControl')
	Component               = kwargs.get('Component', 'atm')
	Bool_Convert_Unit       = kwargs.get('Bool_Convert_Unit', True)
	Bool_LongTermSimulation = kwargs.get('Bool_LongTermSimulation', False)

	Var_kwarg = {\
		'MiPert': MiPert, \
		'Compset': Compset, \
		'Component': Component, \
		'Bool_Convert_Unit': Bool_Convert_Unit, \
		'Bool_LongTermSimulation': Bool_LongTermSimulation, \
	}

	if (Run == 'Gap'):
		
		Suffix_Run = ''
	
	else:
		
		Suffix_Run = '_' + Run

	if ((MiPert is None) or (MiPert == 0)):
		
		Suffix_MiPert = ''

	else:
		
		Suffix_MiPert = '_e' + str(MiPert)

	if (Bool_LongTermSimulation):
		
		if (Compset == 'B2000'):

			File_Path = '/work3/CESM2_DEFExp_piControl/CESM2_DEFExp_{Compset}{Suffix_Run}{Suffix_MiPert}_LongTermSimulation/{Component}/hist/CESM2*.h0.*.nc'.format(Compset=Compset, Suffix_Run=Suffix_Run, Suffix_MiPert=Suffix_MiPert, Component=Component)

		else:

			File_Path = '/work3/CESM2_DEFExp_piControl/CESM2_DEFExp_{Compset}_{Year:04d}{Suffix_Run}{Suffix_MiPert}_LongTermSimulation/{Component}/hist/CESM2*.h0.*.nc'.format(Compset=Compset, Year=Year, Suffix_Run=Suffix_Run, Suffix_MiPert=Suffix_MiPert, Component=Component)

	else:
		
		File_Path = '/work3/CESM2_DEFExp_piControl/CESM2_DEFExp_{Compset}_{Year:04d}{Suffix_Run}{Suffix_MiPert}/{Component}/hist/CESM2*.h0.*.nc'.format(Compset=Compset, Year=Year, Suffix_Run=Suffix_Run, Suffix_MiPert=Suffix_MiPert, Component=Component)

	if ((Run == 'DEF') and (MiPert != 0)):

		return None
	
	if (Var == 'PRECT'):

		Data = Get_RestartSimulation(Run, Year, 'PRECC', **Var_kwarg) + Get_RestartSimulation(Run, Year, 'PRECL', **Var_kwarg)

	elif (Var == 'FLUS'):
		
		# FLUS: FLNS + FLDS
		Data = Get_RestartSimulation(Run, Year, 'FLNS', **Var_kwarg) + Get_RestartSimulation(Run, Year, 'FLDS', **Var_kwarg)

	else:
		
		RawData   = nc.MFDataset(File_Path, 'r')
		Data      = RawData.variables[Var][:]

	if (np.ndim(Data)>=4):

		PS = RawData.variables['PS'][:]
		Data = Calc_VerticalInterpolation(Data, PS=PS)

	if (Bool_Convert_Unit):

		if (Var in ['PRECC', 'PRECL'])  : Data = Data * Lv * 1e3
		elif (Var in ['QSOIL', 'QVEGE', 'QVEGT']): Data = Data * Lv

	return Data

def Get_MSE_Budget(Run, Year, Var, **kwargs):

	MiPert                  = kwargs.get('MiPert', None)
	Compset                 = kwargs.get('Compset', 'piControl')
	Bool_LongTermSimulation = kwargs.get('Bool_LongTermSimulation', False)

	if ((MiPert is None) or (MiPert == 0)):
		
		Suffix_MiPert = ''

	else:
		
		Suffix_MiPert = '_e' + str(MiPert)

	File_Path = 'Output_Data/nc_MSE_Budget/'

	if (Bool_LongTermSimulation):
		
		if (Compset == 'B2000'):

			File_Name = 'nc_MSE_Budget_{Compset}_LongTermSimulation.nc'.format(Compset=Compset)

		else:

			File_Name = 'nc_MSE_Budget_{Compset}_{Year:04d}_LongTermSimulation.nc'.format(Compset=Compset, Year=Year)

	else:
		
		File_Name = 'nc_MSE_Budget_{Compset}_{Year:04d}{Suffix_MiPert}.nc'.format(Compset=Compset, Year=Year, Suffix_MiPert=Suffix_MiPert)
	
	RawData   = nc.Dataset(File_Path + File_Name, 'r')
	Data      = RawData.variables['{}_{}'.format(Var, Run)][:]

	return Data

def Get_Events_List(**kwargs):

	Filter      = kwargs.get('Filter', None)
	File_Dir    = kwargs.get('File_Dir', '/work3/CESM2_DEFExp_piControl/')

	File_List   = [i for i in os.listdir(File_Dir) if ('CESM2_DEFExp_piControl' in i) and (not '_LongTermSimulation' in i)]
	Events_List = []
	
	for i_Event in np.unique(np.array([int(i.split('_')[3]) for i in File_List])):

		if (Filter == 'Entire'):

			if (len([i for i in File_List if '{:04d}'.format(i_Event) in i]) in [12]):

				Events_List.append(i_Event)

		elif (Filter == 'MaPert_Only'):

			if (len([i for i in File_List if '{:04d}'.format(i_Event) in i]) in [3, 4, 12]):

				Events_List.append(i_Event)

		else:

			Events_List.append(i_Event)
	
	return Events_List