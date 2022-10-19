# LinePlot_SeasonalCycle_PRECT_VIMSE.py

# Import Module
import numpy as np
import pandas as pd

import sys

sys.path.append('..')
import ClimPlot.DataProcessing.DataProcessing as DataProc
import ClimPlot.GetData.GetData as GetData
import ClimPlot.Plot.LinePlot_Seasonal as CPLinePlot

def Get_Data_Config(ModelProperties, **kwargs):

	Data_Config = {}
	Data_Config['Year_Start'] = 1970
	Data_Config['Year_End']   = 2004

	Data_Config['n_time']       = (Data_Config['Year_End'] - Data_Config['Year_Start'] + 1) * 12
	Data_Config['n_year']       = (Data_Config['Year_End'] - Data_Config['Year_Start'] + 1)

	Data_Config['Lev_Interp']   = DataProc.Get_Lev_Interp()

	Data_Config['n_Lev_Interp'] = len(Data_Config['Lev_Interp'])

	if ('Ensemble' in ModelProperties):

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
	Table_GetMergeData = Table_GetMergeData.append(pd.DataFrame({'VarName': ['PRECT'], 'VarName_File': ['PRECT'], 'Range': [Main_Data_Range]}))

	Data = {}

	for ind_Var, i_Var in Table_GetMergeData.iterrows():

		Var_Config = {'VarName': i_Var['VarName'], 'VarName_File': i_Var['VarName_File'], 'Range': i_Var['Range']}

		# Set time range
		Var_Config['Data_TimeRange'] = [Data_Config['Year_Start'], Data_Config['Year_End']]
		
		if ('YearShift' in ModelProperties['Properties']):
			
			Var_Config['YearShift'] = ModelProperties['Properties']['YearShift']

		Data[i_Var['VarName']] = GetData.GetMergeData(ModelProperties, Var_Config, Verbose=True)

		for i_Run in ['CTR', 'DEF']:

			# Preprocessing: Calculate spatial mean over land region
			Data[i_Var['VarName']][i_Run] = DataProc.Calc_LonLat_Average(ModelProperties, Data[i_Var['VarName']][i_Run], Range=i_Var['Range'])
			
			if (i_Var['VarName'] in ['MSE']):

				# Vertical integration of MSE
				Integration_Limits = slice(1, -1)
				Data[i_Var['VarName']][i_Run] = -DataProc.Calc_VI(Data[i_Var['VarName']][i_Run][..., Integration_Limits], \
																  np.array(Data_Config['Lev_Interp'])[Integration_Limits] * 100, -1)
				
		Data[i_Var['VarName']]['ANO'] = Data[i_Var['VarName']]['DEF'] - Data[i_Var['VarName']]['CTR']

	# Rename
	Data['VI_MSE'] = Data['MSE']
	del Data['MSE']

	return Data

def Get_YTicks(Plot_Type, VarName, Run):

	YTicks_Dict = {}
	YTicks_Dict['Mean'] = {'PRECT': {'CTR': np.arange(100, 350, 50), 'DEF': np.arange(100, 350, 50), 'ANO': np.arange(-10, 50, 10)}}
	YTicks_Dict['Std']  = {'PRECT': {'CTR': np.arange(0, 100, 20), 'DEF': np.arange(0, 100, 20), 'ANO': np.arange(0, 100, 20)}}
	
	return YTicks_Dict[Plot_Type][VarName][Run]

###########################################

if __name__ == '__main__':
	
	# Set data source
	Model             = 'AMIP'
	Main_Data_Range   = 'MC_data_2'
	ModelProperties   = Get_ModelProperties(Model)
	Data_Config       = Get_Data_Config(ModelProperties)

	# Get data
	Data              = Get_MergeData(ModelProperties)

	for i_Run in ['CTR', 'DEF', 'ANO']:

		for i_VarName, i_Unit in zip(['PRECT'], ['$W\cdot m^{-2}$']):

			Seasonal_Mean, Seasonal_Mean_CI, Seasonal_Variability = DataProc.Calc_SeasonalCycle(Data[i_VarName][i_Run], 1, Return_Variability='Std')
			Seasonal_Mean = np.nanmean(Seasonal_Mean, axis=0)
			
			# Plot: Mean
			YTicks = Get_YTicks('Mean', i_VarName, i_Run)
			Plot_Config = {\
				'VarName'        : '{}_SeasonalCycle_Mean'.format(i_VarName), \
				'VarLongName'    : '{} SeasonalCycle_Mean'.format(i_VarName), \
				'FigureName_Note': '_{}_Mean'.format(i_Run), \
				'YTicks'         : YTicks, \
				'YLim'           : [YTicks[0], YTicks[-1]], \
				'YLabel'         : '{} ({})'.format(i_VarName, i_Unit) \
			}
			
			CPLinePlot.Plot(ModelProperties, {'LinePlot': Seasonal_Mean, 'LinePlot_CI': Seasonal_Mean_CI}, Data_Config=Data_Config, Plot_Config=Plot_Config)

			# Plot: Std
			YTicks = Get_YTicks('Std', i_VarName, i_Run)
			Plot_Config = {\
				'VarName'        : '{}_SeasonalCycle_Std'.format(i_VarName), \
				'VarLongName'    : '{} SeasonalCycle_Std'.format(i_VarName), \
				'FigureName_Note': '_{}_Std'.format(i_Run), \
				'YTicks'         : YTicks, \
				'YLim'           : [YTicks[0], YTicks[-1]], \
				'YLabel'         : '{} Std. ({})'.format(i_VarName, i_Unit) \
			}
			
			CPLinePlot.Plot(ModelProperties, {'LinePlot': Seasonal_Variability}, Data_Config=Data_Config, Plot_Config=Plot_Config)
