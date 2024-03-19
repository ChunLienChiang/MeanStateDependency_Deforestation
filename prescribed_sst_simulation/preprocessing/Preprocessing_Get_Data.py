# Preprocessing_Get_Data.py
import numpy as np
import xarray as xr
import subprocess
import os
import multiprocessing as mp
import functools as ft
import json
import sys
sys.path.append('../')
import preprocessing.Preprocessing as Prep

# Get configuration
Config = json.load(open(os.path.join(os.path.dirname(__file__), '../config.json')))

# Declare constants
Cp = 1006
g  = 9.8
Lv = 2.25e6

def Calc_VerticalInterpolation(Data, PS, hyam, hybm):

	import Ngl

	# Interpolation
	Data_Interp = Ngl.vinth2p(\
		Data, hyam, hybm, \
		Prep.Get_Lev_Interp(), PS, \
		1, 1000, 1, False\
	)
	
	# Mask out invalid values
	Data_Interp = np.where(Data_Interp>1e25, np.nan, Data_Interp)

	return Data_Interp

def _Process_Get_Simulation(File, Var):

	# Connect to nc file
	RawData   = xr.open_dataset(File, engine='netcdf4')
	Data      = RawData[Var].to_numpy()
	
	return Data

def Get_Simulation(Run, En, Var, **kwargs):

	Range             = kwargs.get('Range', 'Global_Analysis')
	Bool_Convert_Unit = kwargs.get('Bool_Convert_Unit', True)
	Component         = kwargs.get('Component', 'atm')
	n_processes       = kwargs.get('n_processes', 8)

	Var_kwarg = {\
		'Component': Component, \
		'Bool_Convert_Unit': Bool_Convert_Unit, \
	}

	# If the argument En is None, set it to En5
	if (En is None):

		En = 5

	# Set file path and file name
	File_Path = Config['Data_Path'] + '{Case}_{En:02d}/{Component}/hist/'.format(Case=Config['Data_Name'][Run], En=En, Component=Component)
	File_List = [\
		File_Path + f for f in os.listdir(File_Path) if \
		(int(f.split('.')[-2].split('-')[0]) >= Config['Data_TimeRange'][0]) and \
		(int(f.split('.')[-2].split('-')[0]) <= Config['Data_TimeRange'][1]) and \
		(f.split('.')[-1] in ['nc', 'nc4'])\
	]
	File_List = sorted(File_List, key=lambda f: f.split('.')[-2])

	if (Var == 'PRECT'):

		# PRECT = PRECC + PRECL
		Data = \
			Get_Simulation(Run, En, 'PRECC', **Var_kwarg) + \
			Get_Simulation(Run, En, 'PRECL', **Var_kwarg)

	elif (Var == 'FLUS'):
		
		# FLUS: FLNS + FLDS
		Data = \
			Get_Simulation(Run, En, 'FLNS', **Var_kwarg) + \
			Get_Simulation(Run, En, 'FLDS', **Var_kwarg)

	else:

		# Create object for multiprocessing
		pool        = mp.Pool(n_processes)
		pool_result = pool.map(ft.partial(_Process_Get_Simulation, Var=Var), File_List)
		pool.close()
		pool.join()
		Data        = np.concatenate(pool_result, axis=0)
	
	# Determine whether the data have vertical dimension and call function for vertical interpolation
	if (np.ndim(Data)>=4):

		PS        = np.concatenate(mp.Pool(n_processes).map(ft.partial(_Process_Get_Simulation, Var='PS'), File_List), axis=0)
		hyam      = Prep.Get_RefData(Var='hyam')
		hybm      = Prep.Get_RefData(Var='hybm')
		Data      = Calc_VerticalInterpolation(Data, PS, hyam, hybm)
	
	if (Bool_Convert_Unit):

		if (Var in ['PRECC', 'PRECL']): Data = Data * Lv * 1e3
		elif (Var in ['QSOIL', 'QVEGE', 'QVEGT']): Data = Data * Lv

	# Crop data
	if (Range != 'Global_Analysis'):

		Data = Prep.Crop_Range(Data, Range)

	return Data

def Get_Simulation_Extract(Run, En, Var, Range, Component='atm'):
	
	File_Name = 'ModelOutput_Extract.{Component}.{Var}.{Run}.{Year_Start}-{Year_End}.En{En}.{Range}.nc'.format(\
		Component=Component, \
		Var=Var, \
		Run=Run, \
		Year_Start=Config['Data_TimeRange'][0], \
		Year_End=Config['Data_TimeRange'][1], \
		En=En, \
		Range=Range, \
	)

	Data = _Process_Get_Simulation(os.path.join(os.path.dirname(__file__), '../output/Output_Data/ModelOutput_Extract/{}'.format(File_Name)), Var)

	return Data

if (__name__ == '__main__'):
	
	"""
	When this file is called, try to get temperature in the case DEF (en-5)
	"""
	
	Data = Get_Simulation('DEF', 5, 'T')
	print('Get temperature data')
	print('Shape:', Data.shape)