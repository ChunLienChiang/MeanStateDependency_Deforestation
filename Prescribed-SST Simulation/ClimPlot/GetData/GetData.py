# GetData.py

# Import modules
import numpy as np
import pandas as pd
import netCDF4 as nc
import os

import sys
sys.path.append('..')
import ClimPlot.DataProcessing.DataProcessing as DataProc

def Process_GetMergeData(i_En, ModelProperties, Verbose, Var_Config, i_Run, Data_Path, Time_Slice, Lat_Slice, Lon_Slice):

	if (Verbose):

		print('Get Data: {}, {}, {}'.format(Var_Config['VarName'], i_Run, str(i_En).zfill(ModelProperties['Ensemble']['Default_En_Format']['ZeroFill'])))

	Data_FileName    = '/{ModelName}_{En}_{VarName}_{Run}.nc'.format(ModelName=ModelProperties['ModelName'], \
																	 En=str(i_En).zfill(ModelProperties['Ensemble']['Default_En_Format']['ZeroFill']), \
																	 VarName=Var_Config['VarName'], Run=i_Run)
	Rawdata          = nc.Dataset(Data_Path + Data_FileName)
	Var_Rawdata      = Rawdata.variables[Var_Config['VarName_File']][Time_Slice, ..., Lat_Slice, Lon_Slice].squeeze()
	
	return Var_Rawdata[None, :]

def Inquire_DataRange(ModelProperties, VarName):
	
	# Inquire data range from Table_VerMerge.csv file
	Table_VarMerge = pd.read_csv('../{ModelPath}/Table_VarMerge.csv'.format(ModelPath=ModelProperties['ModelPath']))
	Origin_Range   = Table_VarMerge[Table_VarMerge['VarName']==VarName]

	if (Origin_Range.empty):

		print('Variable does not exist: {}'.format(VarName))
		quit()

	else:

		return Origin_Range['DataRange'].to_numpy()[0]

def GetMergeData(ModelProperties, Var_Config, **kwargs):

	MultiProcessing     = kwargs.get('MultiProcessing', True)
	MultiProcessing_job = kwargs.get('MultiProcessing_job', 2)
	Verbose             = kwargs.get('Verbose', False)

	# Connect to path
	Data_Path = os.path.abspath('../{ModelPath}/Data/{VarName}/'.format(ModelPath=ModelProperties['ModelPath'], VarName=Var_Config['VarName']))
	
	# Get ensemble list
	if ('En' in Var_Config):

		En = Var_Config['En']
	
	else:

		En = ModelProperties['Ensemble']['Default_En']

	# Calculate spatial range
	if ('Range' in Var_Config):
		
		# Crop spatial range
		Origin_Range = Inquire_DataRange(ModelProperties, Var_Config['VarName'])

		Origin_Lon_Min, Origin_Lon_Max, Origin_Lat_Min, Origin_Lat_Max = DataProc.Get_Range(ModelProperties, Range=Origin_Range, Type='Grid')
		Lon_Min, Lon_Max, Lat_Min, Lat_Max                             = DataProc.Get_Range(ModelProperties, Range=Var_Config['Range'], Type='Grid')

		Lon_Min   = Lon_Min - Origin_Lon_Min
		Lon_Max   = Lon_Max - Origin_Lon_Max
		Lat_Min   = Lat_Min - Origin_Lat_Min
		Lat_Max   = Lat_Max - Origin_Lat_Max
		
		Lon_Max = None if Lon_Max == 0 else Lon_Max
		Lat_Max = None if Lat_Max == 0 else Lat_Max

		Lat_Slice = slice(Lat_Min, Lat_Max, 1)
		Lon_Slice = slice(Lon_Min, Lon_Max, 1)

	else:

		Lat_Slice = slice(None)
		Lon_Slice = slice(None)
	
	# Calculate temporal range
	if ('Data_TimeRange' in Var_Config):
		
		if (Var_Config['Data_TimeRange'][0]!=ModelProperties['Properties']['Data_TimeRange'][0]) or \
		   (Var_Config['Data_TimeRange'][-1]!=ModelProperties['Properties']['Data_TimeRange'][1]):
			
			Data_TimeRange = [Var_Config['Data_TimeRange'][0] - ModelProperties['Properties']['Data_TimeRange'][0], \
							  Var_Config['Data_TimeRange'][-1] - ModelProperties['Properties']['Data_TimeRange'][0]]
			Time_Slice     = slice((Data_TimeRange[0])*12, (Data_TimeRange[-1]+1)*12, 1)
		
		else:

			Time_Slice     = slice(None)

	else:
		
		Time_Slice     = slice(None)
	
	# Get the run which have to be obtained
	if ('ModelRun' in Var_Config):

		Run = Var_Config['ModelRun']

	else:

		Run = ModelProperties['ModelRun']

	# Create new dictionary to store data
	Data = {}

	for i_Run in Run:

		Data[i_Run] = []

		if (MultiProcessing):
			
			import multiprocessing as mp
			import functools as fts
			
			Process_GetMergeData_Partial = fts.partial(Process_GetMergeData, ModelProperties=ModelProperties, Verbose=Verbose, Var_Config=Var_Config, i_Run=i_Run, \
													   Data_Path=Data_Path, Time_Slice=Time_Slice, Lat_Slice=Lat_Slice, Lon_Slice=Lon_Slice)

			MP_pool     = mp.Pool(processes=MultiProcessing_job)
			MP_Result   = MP_pool.map(Process_GetMergeData_Partial, En)
			
			Data[i_Run] = np.concatenate(tuple(MP_Result), axis=0)
		
		else:
			
			for i_En in En:

				print('Get Data: {}, {}, {}'.format(Var_Config['VarName'], i_Run, str(i_En).zfill(ModelProperties['Ensemble']['Default_En_Format']['ZeroFill'])))

				Data_FileName    = '/{ModelName}_{En}_{VarName}_{Run}.nc'.format(ModelName=ModelProperties['ModelName'], \
																				En=str(i_En).zfill(ModelProperties['Ensemble']['Default_En_Format']['ZeroFill']), \
																				VarName=Var_Config['VarName'], Run=i_Run)
				Rawdata          = nc.Dataset(Data_Path + Data_FileName)
				Var_Rawdata      = Rawdata.variables[Var_Config['VarName_File']][Time_Slice, ..., Lat_Slice, Lon_Slice].squeeze()
				Data[i_Run].append(Var_Rawdata[None, :])

			Data[i_Run] = np.concatenate(tuple(Data[i_Run]), axis=0)

	return Data