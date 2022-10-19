# DataAnalysis_PRECT.py

# Import Module
import numpy as np
import pandas as pd
import os

import sys
sys.path.append('..')
import ClimPlot.DataProcessing as DataProc
import ClimPlot.GetData.GetData as GetData
import ClimPlot.Plot.ScatterPlot as CPScatter

def Get_Data_Config(ModelProperties, **kwargs):

	Data_Config = {}
	Data_Config['Year_Start'] = 1970
	Data_Config['Year_End']   = 2004

	Data_Config['n_time']       = (Data_Config['Year_End'] - Data_Config['Year_Start'] + 1) * 12
	Data_Config['n_year']       = (Data_Config['Year_End'] - Data_Config['Year_Start'] + 1)

	if ('Ensemble' in ModelProperties):

		Data_Config['n_En']       = len(ModelProperties['Ensemble']['Default_En'])
		Data_Config['Dim_Time']   = 1

	else:

		Data_Config['Dim_Time']   = 0

	return Data_Config

def Get_ModelProperties(Model):

	import json
	Json_File            = open('ModelProperties.json')
	ModelProperties      = json.load(Json_File)['Model'][Model]

	return ModelProperties

def Get_MergeData(ModelProperties, **kwargs):

	Table_GetMergeData = pd.DataFrame({'VarName': [], 'VarName_File': [], 'Range': []})
	Table_GetMergeData = Table_GetMergeData.append(pd.DataFrame({'VarName': ['PRECT'], 'VarName_File': ['PRECT'], 'Range': [Main_Data_Range]}))
	Table_GetMergeData = Table_GetMergeData.append(pd.DataFrame({'VarName': ['PRECC'], 'VarName_File': ['PRECC'], 'Range': [Main_Data_Range]}))

	Data = {}

	for ind_Var, i_Var in Table_GetMergeData.iterrows():

		Var_Config = {'VarName': i_Var['VarName'], 'VarName_File': i_Var['VarName_File'], 'Range': i_Var['Range']}

		# Set time range
		Var_Config['Data_TimeRange'] = [Data_Config['Year_Start'], Data_Config['Year_End']]
		
		if ('YearShift' in ModelProperties['Properties']):
			
			Var_Config['YearShift'] = ModelProperties['Properties']['YearShift']

		Data[i_Var['VarName']] = GetData.GetMergeData(ModelProperties, Var_Config)

		for i_Run in ['CTR', 'DEF']:

			# Preprocessing: Calculate spatial mean over land region
			Data[i_Var['VarName']][i_Run] = DataProc.DataProcessing.Calc_LonLat_Average(ModelProperties, Data[i_Var['VarName']][i_Run], Range=i_Var['Range'])
			
			# Preprocessing: Mask data by month range
			Data[i_Var['VarName']][i_Run] = DataProc.DataProcessing.Calc_Month_Mask(Data[i_Var['VarName']][i_Run], i_Month, Data_Config['Dim_Time'])

			# Preprocessing: Calculate annual mean
			Data[i_Var['VarName']][i_Run] = DataProc.DataProcessing.Calc_AnnualMean(Data[i_Var['VarName']][i_Run], Data_Config['Dim_Time'])
		
		Data[i_Var['VarName']]['ANO'] = Data[i_Var['VarName']]['DEF'] - Data[i_Var['VarName']]['CTR']

	return Data

def Get_Grouping(Month, **kwargs):

	GroupBy       = kwargs.get('GroupBy', 'Percentile')

	if (GroupBy=='Percentile'):

		Grouping = pd.read_csv('Output_Result/{}/MSEProfile_Group/DataFrame_PercentileRange_CTR.csv'.format(ModelProperties['ModelName']))
		Grouping = Grouping[Grouping['Month']==Month]
		Grouping = Grouping[['Year', 'En', 'Percentile_Group', 'ANO_P']].sort_values(by=['En', 'Year'])
		Grouping = Grouping['Percentile_Group'].to_numpy().reshape((-1, Data_Config['n_year']))
	
	return Grouping

def SaveFigure(Data):
	
	print('Save Figures')
	
	for i_Month in Month_List:

		i_Data         = Data[Data['Month']==i_Month]

		# Scatter plot
	
		Plot_Config = {\
			'VarName'         : 'PRECT_CTR_ANO', \
			'VarLongName'     : 'Precipitation CTR~ANO', \
			'FigureName_Note' : '_{}'.format(i_Month), \
			'XTicks'          : np.arange(0, 400, 100), \
			'XLabel'          : r'CTL PRECT ($W\cdot m^{-2}$)', \
			'YTicks'          : np.arange(-50, 150, 50), \
			'YLabel'          : r'ANO PRECT ($W\cdot m^{-2}$)', \
			'LinRegLine'      : True \
		}

		CPScatter.Plot(ModelProperties, {'X': i_Data['PRECT_CTR'], 'Y': i_Data['PRECT_ANO']}, \
					   Data_Config=Data_Config, Plot_Config=Plot_Config)

	return
	
###########################################

if __name__ == '__main__':

	import scipy

	# Set data source
	Model             = 'AMIP'
	Main_Data_Range   = 'MC_data_2'
	ModelProperties   = Get_ModelProperties(Model)
	Data_Config       = Get_Data_Config(ModelProperties)
	Output_Reg        = False

	DataFrame = pd.DataFrame({'Month': [], 'PRECT_CTR': [], 'PRECT_DEF': [], 'PRECT_ANO': [], 'Group': []})
	
	Month_List        = ['All', 'NDJ', 'JA', 'DJF', 'January', 'July']

	for i_Month in Month_List:
	
		# Get data
		Data              = Get_MergeData(ModelProperties)
		Grouping          = Get_Grouping(i_Month)

		DataFrame = DataFrame.append(pd.DataFrame({'Month': [i_Month] * Data_Config['n_En'] * Data_Config['n_year'], \
												   'PRECT_CTR': Data['PRECT']['CTR'].ravel(), \
												   'PRECT_DEF': Data['PRECT']['DEF'].ravel(), \
												   'PRECT_ANO': Data['PRECT']['ANO'].ravel(), \
												   'Group': Grouping.ravel()}))
		
		Mask_10 = (Data['PRECT']['CTR'].ravel()<=np.nanpercentile(Data['PRECT']['CTR'].ravel(), 10))
		Mask_90 = (Data['PRECT']['CTR'].ravel()>=np.nanpercentile(Data['PRECT']['CTR'].ravel(), 90))
		print(i_Month)
		print(scipy.stats.percentileofscore(Data['PRECT']['DEF'].ravel(), np.nanmean(Data['PRECT']['DEF'].ravel()[Mask_10])))
		print(scipy.stats.percentileofscore(Data['PRECT']['DEF'].ravel(), np.nanmean(Data['PRECT']['DEF'].ravel()[Mask_90])))

	# Save figures
	SaveFigure(DataFrame)

	if not os.path.exists('Output_Result/{}/DataAnalysis_PRECT/'.format(ModelProperties['ModelName'])):

		os.mkdir('Output_Result/{}/DataAnalysis_PRECT/'.format(ModelProperties['ModelName']))
	
	DataFrame.to_csv('Output_Result/{}/DataAnalysis_PRECT/DataAnalysis_PRECT.csv'.format(ModelProperties['ModelName']), index=False)