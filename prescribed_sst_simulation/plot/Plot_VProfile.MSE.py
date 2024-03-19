"""
Plot_VProfile.MSE.py
==========================
Plot vertical profiles for the MSE.
"""

import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.ticker as mplticker
import scipy.stats as scista
import json
import os
import sys
sys.path.append('../')
import preprocessing.Preprocessing as Prep

# Get configuration
Config = json.load(open(os.path.join(os.path.dirname(__file__), '../config.json')))

def Plot_VProfile(Plot_Data):

	# Create matplotlib object
	fig, ax = plt.subplots(figsize=(4, 4), dpi=300)

	# Plot: vertical profiles of MSE
	ax.plot(\
		Plot_Data['Mean_Group1'], Prep.Get_Lev_Interp(), \
		linewidth=2.5, color='C0', \
		label='{PrecipMean:.1f} '.format(PrecipMean=Plot_Data['PrecipMean_Group1']) + r'$W\cdot m^{-2}$' + '\n~1st quartile\nn={n:2d}'.format(n=Plot_Data['n_Group1']), \
	)
	ax.plot(\
		Plot_Data['Mean_Group2'], Prep.Get_Lev_Interp(), \
		linewidth=2.5, color='C3', \
		label='{PrecipMean:.1f} '.format(PrecipMean=Plot_Data['PrecipMean_Group2']) + r'$W\cdot m^{-2}$' + '\n3th quartile~\nn={n:2d}'.format(n=Plot_Data['n_Group2']), \
	)

	# Plot: confidence interval
	ax.fill_betweenx(\
		Prep.Get_Lev_Interp(), Plot_Data['CI_Group1'][0], Plot_Data['CI_Group1'][1], \
		color='C0', ec='none', alpha=0.15, \
	)
	ax.fill_betweenx(\
		Prep.Get_Lev_Interp(), Plot_Data['CI_Group2'][0], Plot_Data['CI_Group2'][1], \
		color='C3', ec='none', alpha=0.15, \
	)

	# Plot: symbols for significance test
	ax.plot(\
		np.array([-390]*21), \
		np.where((Plot_Data['pvalues_Diff']<=0.10)&(Plot_Data['pvalues_Diff']>0.05), Prep.Get_Lev_Interp(), np.nan), \
		linestyle='None', \
		marker='o', \
		markersize=5, \
		markerfacecolor='Grey', \
		markeredgecolor='none', \
	)
	ax.plot(\
		np.array([-390]*21), \
		np.where(Plot_Data['pvalues_Diff']<0.05, Prep.Get_Lev_Interp(), np.nan), \
		linestyle='None', \
		marker='<', \
		markersize=5, \
		markerfacecolor='Grey', \
		markeredgecolor='none', \
	)

	# Plot: reference line
	ax.vlines(0, 0, 1000, colors='Grey', linewidth=1.5, zorder=-2)

	# Configuration
	ax.set_xlabel(r'MSE ($J\cdot kg^{-1}$)')
	ax.set_xlim([-400, 400])
	ax.set_xticks(np.arange(-400, 600, 200))
	ax.set_ylabel('Pressure (hPa)')
	ax.set_ylim([975, 400])
	ax.set_yticks([975, 900, 800, 700, 600, 500, 400])
	ax.tick_params(axis='y')

	ax.yaxis.set_major_formatter(mplticker.ScalarFormatter())
	ax.yaxis.get_major_formatter().set_scientific(False)
	ax.yaxis.get_major_formatter().set_useOffset(False)

	# Legend
	ax.legend(loc='upper right', fontsize=7)
	
	# Plot save
	plt.tight_layout()
	Figure_Path = '../output/Output_Figure/Plot_VProfile.MSE/'
	Figure_Name = 'Plot_VProfile.MSE.png'
	if (not os.path.exists(Figure_Path)): os.makedirs(Figure_Path)
	plt.savefig(Figure_Path + Figure_Name)
	plt.close('all')

	return

if (__name__ == '__main__'):

	# Read data: grouping
	Data_PRECT  = pd.read_csv('../output/Output_Data/Group/Group.1970-2005.csv')
	Mask_Group1 = (Data_PRECT['Group']==0)
	Mask_Group2 = (Data_PRECT['Group']==3)

	# Read data: MSE
	Data_MSE = []

	for i_En in Config['Data_En']:

		Data_File_Path = '../output/Output_Data/MSE_Budget/'
		Data_File_Name = 'MSE_Budget.1970-2005.En{En}.MC.nc'.format(En=i_En)
		if not (os.path.exists(Data_File_Path + Data_File_Name)): raise ValueError('The file {} does not exist.'.format(Data_File_Name))

		Data_MSE.append(xr.open_dataset(Data_File_Path + Data_File_Name)['MSE_CTL'].to_numpy())

	Data_MSE = np.concatenate(tuple(Data_MSE), axis=0)
	
	# Calculate anomalies of MSE vertical profile
	Data_MSE = Data_MSE - np.nanmean(Data_MSE, axis=0)

	# Plot: VProfile
	Plot_Data = {\
		'Mean_Group1'      : Prep.Calc_CI(Data_MSE[Mask_Group1, :], Data_Type='General', Test_Axis=0)[0], \
		'CI_Group1'        : Prep.Calc_CI(Data_MSE[Mask_Group1, :], Data_Type='General', Test_Axis=0)[1], \
		'Mean_Group2'      : Prep.Calc_CI(Data_MSE[Mask_Group2, :], Data_Type='General', Test_Axis=0)[0], \
		'CI_Group2'        : Prep.Calc_CI(Data_MSE[Mask_Group2, :], Data_Type='General', Test_Axis=0)[1], \
		'PrecipMean_Group1': np.nanmean(Data_PRECT.loc[Mask_Group1, 'PRECT_ANO'].to_numpy()), \
		'PrecipMean_Group2': np.nanmean(Data_PRECT.loc[Mask_Group2, 'PRECT_ANO'].to_numpy()), \
		'n_Group1'         : np.sum(Mask_Group1), \
		'n_Group2'         : np.sum(Mask_Group2), \
		'pvalues_Diff'     : scista.ttest_ind(Data_MSE[Mask_Group1, :], Data_MSE[Mask_Group2, :], axis=0, equal_var=False, nan_policy='omit', alternative='greater').pvalue
	}

	Plot_VProfile(Plot_Data)