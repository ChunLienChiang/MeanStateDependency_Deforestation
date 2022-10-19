# VProfile_MSEProfile_PercentileRange.py

# Import Module
import numpy as np
import pandas as pd
import scipy.stats as SciSt
import sys

sys.path.append('..')
import ClimPlot.DataProcessing as DataProc
import ClimPlot.GetData.GetData as GetData
import ClimPlot.Plot.VProfile as CPVProfile

def Get_Data_Config(ModelProperties, **kwargs):

	Data_Config = {}
	Data_Config['Year_Start'] = ModelProperties['Properties']['Data_TimeRange'][0]
	Data_Config['Year_End']   = ModelProperties['Properties']['Data_TimeRange'][-1]

	Data_Config['n_time']       = (Data_Config['Year_End'] - Data_Config['Year_Start'] + 1) * 12
	Data_Config['n_year']       = (Data_Config['Year_End'] - Data_Config['Year_Start'] + 1)

	Data_Config['Lev_Interp']   = DataProc.DataProcessing.Get_Lev_Interp()

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
	Table_GetMergeData = pd.concat( \
		[Table_GetMergeData, pd.DataFrame({'VarName': ['MSE'], 'VarName_File': ['MSE'], 'Range': [Main_Data_Range]}), pd.DataFrame({'VarName': ['PRECT'], 'VarName_File': ['PRECT'], 'Range': [Main_Data_Range]})], \
		axis=0
	)

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
			
		Data[i_Var['VarName']]['ANO'] = Data[i_Var['VarName']]['DEF'] - Data[i_Var['VarName']]['CTR']

	for ind_Var, i_Var in Table_GetMergeData.iterrows():

		for i_Run in ['CTR', 'DEF', 'ANO']:

			if (i_Var['VarName']=='MSE'):

				# Calculate Calculate climatology anomaly for MSE
				Data[i_Var['VarName']][i_Run] = DataProc.DataProcessing.Calc_ClimAnomaly(Data[i_Var['VarName']][i_Run], Data_Config['Dim_Time'])
			
			# Preprocessing: Mask data by month range
			Data[i_Var['VarName']][i_Run] = DataProc.DataProcessing.Calc_Month_Mask(Data[i_Var['VarName']][i_Run], i_Month, Data_Config['Dim_Time'])
			
			# Preprocessing: Calculate annual mean
			Data[i_Var['VarName']][i_Run] = DataProc.DataProcessing.Calc_AnnualMean(Data[i_Var['VarName']][i_Run], Data_Config['Dim_Time'])

	return Data

###########################################

if __name__ == '__main__':

	Model             = 'AMIP'
	Main_Data_Range   = 'MC_data_2'
	ModelProperties   = Get_ModelProperties(Model)
	Data_Config       = Get_Data_Config(ModelProperties)

	for i_Month in ['All', 'NDJ', 'JA', 'DJF', 'January', 'July']:
		
		# Get data
		Data              = Get_MergeData(ModelProperties)

		# ################################################################

		# Range list
		n_Group         = 4
		Range_List      = np.array([0, 100/n_Group, 100-100/n_Group, 100])

		for i_Run in ['CTR', 'DEF', 'ANO']:

			Range_DataFrame   = pd.concat([pd.DataFrame(Data['PRECT']['ANO'].ravel(), columns=['ANO_P']), \
										   pd.DataFrame(Data['MSE'][i_Run].reshape(-1, Data_Config['n_Lev_Interp']), \
														columns=[str(i_Run)+'_MSE_'+str(i) for i in range(Data_Config['n_Lev_Interp'])])], \
										  axis=1)
			
			# Calculate grouped mean, std, and CI
			Range_DataFrame_Group_Mean  = []
			Range_DataFrame_Group_STD   = []
			Range_DataFrame_Group_n     = []
			Range_DataFrame_Group_CI    = []
			Range_DataFrame_Group_PRECT = []

			for i_Range in range(len(Range_List)-1):

				Range_Min, Range_Max    = np.round(Range_List[i_Range], 1), np.round(Range_List[i_Range+1], 1)
				Range_PRECT_Min, \
				Range_PRECT_Max         = np.percentile(Range_DataFrame['ANO_P'], Range_Min), np.percentile(Range_DataFrame['ANO_P'], Range_Max)

				i_Range_DataFrame       = Range_DataFrame[(Range_DataFrame['ANO_P']>=Range_PRECT_Min)&(Range_DataFrame['ANO_P']<Range_PRECT_Max)].filter(regex='_MSE_')
				i_Range_DataFrame_PRECT = Range_DataFrame[(Range_DataFrame['ANO_P']>=Range_PRECT_Min)&(Range_DataFrame['ANO_P']<Range_PRECT_Max)].mean()['ANO_P']
				i_Range_DataFrame_Mean  = i_Range_DataFrame.mean(axis=0)
				i_Range_DataFrame_STD   = i_Range_DataFrame.std(axis=0)
				i_Range_DataFrame_n     = i_Range_DataFrame.shape[0]

				i_Range_DataFrame_CI    = [(i_Range_DataFrame_Mean - 1.96 * i_Range_DataFrame_STD / np.sqrt(i_Range_DataFrame_n)).to_numpy(), \
										   (i_Range_DataFrame_Mean + 1.96 * i_Range_DataFrame_STD / np.sqrt(i_Range_DataFrame_n)).to_numpy()]

				Range_DataFrame_Group_Mean.append(i_Range_DataFrame_Mean)
				Range_DataFrame_Group_STD.append(i_Range_DataFrame_STD)
				Range_DataFrame_Group_n.append(i_Range_DataFrame_n)
				Range_DataFrame_Group_CI.append(i_Range_DataFrame_CI)
				Range_DataFrame_Group_PRECT.append(i_Range_DataFrame_PRECT)

				Range_DataFrame_Group_New = {'Data': [i_Range_DataFrame_Mean.to_numpy()], 'CI': [i_Range_DataFrame_CI]}

			# Calculate significance level
			Diff_Pvalue = SciSt.ttest_ind_from_stats(\
				Range_DataFrame_Group_Mean[0], Range_DataFrame_Group_STD[0], Range_DataFrame_Group_n[0], \
				Range_DataFrame_Group_Mean[-1], Range_DataFrame_Group_STD[-1], Range_DataFrame_Group_n[-1], equal_var=False, alternative='greater', \
			).pvalue
			SignificanceLevel_095 = np.where(Diff_Pvalue<0.05, True, False)
			SignificanceLevel_090 = np.where(Diff_Pvalue<0.10, True, False)

			# Plot
			Plot_Config = {\
				'VarName'        : 'MSE_Profile_PercentileRange', \
				'VarLongName'    : 'MSE Profile Percentile Range', \
				'FigureName_Note': '_{}_{}_All'.format(i_Run, i_Month), \
				'XTicks'         : np.arange(-400, 600, 200), \
				'XLabel'         : 'MSE', \
				'XUnit'          : r'$J\cdot kg^{-1}$', \
				'Legend_Labels'  : [str(np.round(Range_DataFrame_Group_PRECT[i], 1)) + r' $W\cdot m^{-2}$' + ' (pct.{}~{}), n={}'.format(int(Range_List[i]), int(Range_List[i+1]), Range_DataFrame_Group_n[i]) \
									for i in range(len(Range_List)-1)] \
			}
			
			Plot_Data = {\
				'Data': Range_DataFrame_Group_Mean, \
				'Data_CI': Range_DataFrame_Group_CI, \
				'Data_Color': ['C0', 'Grey', 'C3'], \
				'Data_Sig': {'<': np.where(SignificanceLevel_095, -390, np.nan), 'o': np.where(SignificanceLevel_090&(~SignificanceLevel_095), -390, np.nan)}, \
			}
			CPVProfile.Plot(ModelProperties, Plot_Data, Data_Config=Data_Config, Plot_Config=Plot_Config)