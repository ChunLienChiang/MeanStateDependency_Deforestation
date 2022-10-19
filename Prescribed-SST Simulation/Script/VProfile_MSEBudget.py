# VProfile_MSEBudget.py

# Import Module
import numpy as np
import pandas as pd
import warnings
import sys

sys.path.append('..')
import ClimPlot.DataProcessing.DataProcessing as DataProc
import ClimPlot.Plot.VProfile as CPVProfile

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

def Get_Data(ModelProperties):

	import netCDF4 as nc

	Data = {}

	for i_VarName in ['Omega', 'NOmega', 'hm', 'Vdhm', 'Omegadhm', 'NOmegadhm', 'dV', 'dOmega', 'Vdsd', 'Omegadsd', 'VdQ', 'OmegadQ']:

		Data[i_VarName] = {}
		
		for i_Run in ['CTR', 'DEF', 'ANO', 'ANO_hm', 'ANO_sd', 'ANO_V', 'ANO_Omega', 'ANO_T', 'ANO_Q', 'ANO_NonLinear']:

			# Open File
			Rawdata = nc.Dataset('Output_Result/{}/DataAnalysis_MSEBudget/MSEBudget_{}.nc'.format(ModelProperties['ModelName'], i_Run))
			
			if (i_VarName in Rawdata.variables):
				
				Data[i_VarName][i_Run] = Rawdata.variables[i_VarName][:]
				
				Data[i_VarName][i_Run] = DataProc.Calc_LonLat_Average(ModelProperties, Data[i_VarName][i_Run], Range=Main_Data_Range, Area='Land')
				
				# Preprocessing: Mask data by month range
				Data[i_VarName][i_Run] = DataProc.Calc_Month_Mask(Data[i_VarName][i_Run], Month_Range, Data_Config['Dim_Time'])
				
				# Preprocessing: Calculate annual mean
				Data[i_VarName][i_Run] = DataProc.Calc_AnnualMean(Data[i_VarName][i_Run], Data_Config['Dim_Time'])

	return Data

def Get_XTicks(Plot_Type):

	XTicks = {}

	if (Plot_Type=='Mean'):

		XTicks['Omega']     = {'CTR': np.arange(-5e-2, 7.5e-2, 2.5e-2), 'DEF': np.arange(-5e-2, 7.5e-2, 2.5e-2), 'ANO': np.arange(-2e-2, 3e-2, 1e-2)}
		XTicks['NOmega']    = {'CTR': np.arange(-0.2, 1.2, 0.2), 'DEF': np.arange(-0.2, 1.2, 0.2), 'ANO': np.arange(-0.3, 0.4, 0.1)}
		XTicks['hm']        = {'CTR': np.arange(320000, 370000, 10000), 'DEF': np.arange(320000, 370000, 10000), 'ANO': np.arange(-2000, 3000, 1000)}
		XTicks['Vdhm']      = {'CTR': np.arange(-1e-2, 1.5e-2, 5e-3), 'DEF': np.arange(-1e-2, 1.5e-2, 5e-3), 'ANO': np.arange(-1e-2, 1.5e-2, 5e-3)}
		XTicks['Omegadhm']  = {'CTR': np.arange(-2e-2, 3e-2, 1e-2), 'DEF': np.arange(-2e-2, 3e-2, 1e-2), 'ANO': np.arange(-1e-2, 1.5e-2, 5e-3)}
		XTicks['NOmegadhm'] = {'CTR': np.arange(-2e-2, 3e-2, 1e-2), 'DEF': np.arange(-2e-2, 3e-2, 1e-2), 'ANO': np.arange(-1e-2, 1.5e-2, 5e-3)}
		XTicks['dV']        = {'CTR': np.arange(-1e-5, 1.5e-5, 5e-6), 'DEF': np.arange(-1e-5, 1.5e-5, 5e-6), 'ANO': np.arange(-5e-6, 7.5e-6, 2.5e-6)}
		XTicks['dOmega']    = {'CTR': np.arange(-1e-5, 1.5e-5, 5e-6), 'DEF': np.arange(-1e-5, 1.5e-5, 5e-6), 'ANO': np.arange(-5e-6, 7.5e-6, 2.5e-6)}
		XTicks['Vdsd']      = {'CTR': np.arange(-1e-2, 1.5e-2, 5e-3), 'DEF': np.arange(-1e-2, 1.5e-2, 5e-3), 'ANO': np.arange(-5e-3, 7.5e-3, 2.5e-3)}
		XTicks['Omegadsd']  = {'CTR': np.arange(-5e-2, 7.5e-2, 2.5e-2), 'DEF': np.arange(-5e-2, 7.5e-2, 2.5e-2), 'ANO': np.arange(-2e-2, 3e-2, 1e-2)}
		XTicks['VdQ']       = {'CTR': np.arange(-1e-2, 1.5e-2, 5e-3), 'DEF': np.arange(-1e-2, 1.5e-2, 5e-3), 'ANO': np.arange(-5e-3, 7.5e-3, 2.5e-3)}
		XTicks['OmegadQ']   = {'CTR': np.arange(-5e-2, 7.5e-2, 2.5e-2), 'DEF': np.arange(-5e-2, 7.5e-2, 2.5e-2), 'ANO': np.arange(-2e-2, 3e-2, 1e-2)}
	
	return XTicks

def Get_Grouping(**kwargs):

	GroupBy       = kwargs.get('GroupBy', 'Percentile')

	if (GroupBy=='Percentile'):

		Grouping = pd.read_csv('Output_Result/{}/MSEProfile_Group/DataFrame_PercentileRange_CTR.csv'.format(ModelProperties['ModelName']))
		Grouping = Grouping[Grouping['Month']==Month_Range]
		Grouping = Grouping[['Year', 'En', 'Percentile_Group']].sort_values(by=['En', 'Year'])
		Grouping = Grouping['Percentile_Group'].to_numpy().reshape((-1, Data_Config['n_year']))
	
	return Grouping

###########################################

if __name__ == '__main__':

	warnings.simplefilter('ignore', category=RuntimeWarning)
	
	# Set data source
	Model             = 'AMIP'
	Main_Data_Range   = 'MC_data_2'
	ModelProperties   = Get_ModelProperties(Model)
	Data_Config       = Get_Data_Config(ModelProperties)

	for i_Month in ['All', 'NDJ', 'JA', 'DJF', 'January', 'July']:
	
		# Get data
		Month_Range       = i_Month
		Data              = Get_Data(ModelProperties)
		Grouping          = Get_Grouping()
		
		for i_VarName in ['Omega', 'NOmega', 'hm', 'Vdhm', 'Omegadhm', 'NOmegadhm', 'dV', 'dOmega', 'Vdsd', 'Omegadsd', 'VdQ', 'OmegadQ']:

			for i_Run in ['CTR', 'DEF', 'ANO', 'ANO_hm', 'ANO_sd', 'ANO_V', 'ANO_Omega', 'ANO_T', 'ANO_Q', 'ANO_NonLinear']:
				
				if not (i_Run in Data[i_VarName]):

					continue
				
				if (i_VarName in ['Vdhm', 'Omegadhm', 'NOmegadhm', 'dV', 'dOmega', 'Vdsd', 'Omegadsd', 'VdQ', 'OmegadQ']):

					Data[i_VarName][i_Run] = -Data[i_VarName][i_Run]
				
				Plot_Mean, Plot_CI = [], []

				#for i_Group in np.unique(Grouping):
				for i_Group in [[0], [1, 2], [3]]:

					Group_Mean, Group_CI = DataProc.Calc_CI(np.where(np.isin(Grouping, i_Group)[..., None], Data[i_VarName][i_Run], np.nan), Calc_Type='General')
					Plot_Mean.append(Group_Mean)
					Plot_CI.append(Group_CI)
				
				# Plot: Mean
				Legend_Labels     = [' '] * len(Plot_Mean)
				Legend_Labels[0]  = 'Smallest ANO P'
				Legend_Labels[-1] = 'Greatest ANO P'

				Plot_Config = {\
					'VarName'        : 'MSEBudget_{}'.format(i_VarName), \
					'VarLongName'    : 'MSEBudget {}'.format(i_VarName), \
					'FigureName_Note': '_{}_{}'.format(i_Run, Month_Range), \
					'Legend_Labels'  : Legend_Labels, \
					'XTicks'         : Get_XTicks('Mean')[i_VarName][i_Run[:3]], \
					'XLabel'         : i_VarName \
				}
				
				CPVProfile.Plot(ModelProperties, {'Data': Plot_Mean, 'Data_CI': Plot_CI, 'Data_Color': ['C0', 'Grey', 'C3']}, \
								Data_Config=Data_Config, Plot_Config=Plot_Config)
