"""
Plot_Scatterplot.GMS_Precip.py
==========================
Plot scatterplot of GMS and precipitation
"""

import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import statsmodels.api as sm
import json
import os

# Get configuration
Config = json.load(open(os.path.join(os.path.dirname(__file__), '../config.json')))

def Plot_Scatterplot(Plot_Data, Plot_Config):

	# Create matplotlib object
	fig, ax = plt.subplots(figsize=(4, 4), dpi=300)

	# Plot: scatterplot
	ax.scatter(Plot_Data['Data_Precip'], Plot_Data['Data_GMS'], s=3, color='Blue', alpha=0.9)

	# Plot: regression line
	ax.plot(\
		np.array(Plot_Config['Plot_xlim']), \
		Plot_Data['Linreg_Coef'][0] + Plot_Data['Linreg_Coef'][1]*np.array(Plot_Config['Plot_xlim']), \
		linewidth=1.5, \
		color='Black', \
		alpha=0.9, \
	)

	# Plot: text illustrating the regression line, the r-squared, the slope coefficient and its p-value
	ax.text(\
		0.98, 0.98, \
		r'R$^2$={r2:.3f}'.format(r2=Plot_Data['Linreg_R2']) + '\n' + \
		'p-value={p:.4f}'.format(p=Plot_Data['Linreg_pvalue'][1]) + '\n' + \
		'slope={slope:.3f}'.format(slope=Plot_Data['Linreg_Coef'][1]) + '\n', \
		transform=ax.transAxes, \
		verticalalignment='top', \
		horizontalalignment='right', \
	)

	# Configuration
	ax.set_xlabel('{Run} Precipitation '.format(Run=Plot_Config['Run']) + r'($W\cdot m^{-2}$)')
	ax.set_xlim(Plot_Config['Plot_xlim'])
	ax.set_xticks(np.linspace(*Plot_Config['Plot_xlim'], 4))
	ax.set_ylabel('{Run} GMS '.format(Run=Plot_Config['Run']) + r'($J\cdot kg^{-1}\cdot Pa^{-1}$)')
	ax.set_ylim(Plot_Config['Plot_ylim'])
	ax.set_yticks(np.linspace(*Plot_Config['Plot_ylim'], 5))

	# Plot save
	plt.tight_layout()
	Figure_Path = '../output/Output_Figure/Plot_Scatterplot.GMS_Precip/'
	Figure_Name = 'Plot_Scatterplot.GMS_Precip.{Run}.png'.format(Run=Plot_Config['Run'])
	if (not os.path.exists(Figure_Path)): os.makedirs(Figure_Path)
	plt.savefig(Figure_Path + Figure_Name)
	plt.close('all')

	return

def Get_xlim(Run):
	
	if (Run == 'CTL'):

		xlim = [125, 275]
	
	elif (Run == 'DEF'):

		xlim = [125, 275]
	
	elif (Run == 'ANO'):

		xlim = [-30, 60]
	
	return xlim

def Get_ylim(Run):
	
	if (Run == 'CTL'):

		ylim = [-0.16, 0.16]
	
	elif (Run == 'DEF'):

		ylim = [-0.16, 0.16]
	
	elif (Run == 'ANO'):

		ylim = [-0.08, 0.08]
	
	return ylim

if (__name__ == '__main__'):

	for i_Run in ['CTL', 'DEF', 'ANO']:

		# Get data: precipitation
		Data_Precip = pd.read_csv('../output/Output_Data/Group/Group.1970-2005.csv')['PRECT_{Run}'.format(Run=i_Run)].to_numpy()

		# Get data: GMS
		Data_GMS = np.concatenate([xr.open_dataset('../output/Output_Data/GMS/GMS.1970-2005.En{En}.MC_Extended.nc'.format(En=i), engine='netcdf4')['GMS_{Run}'.format(Run=i_Run)] for i in Config['Data_En']])
		
		# Calculate regression line
		Linreg        = sm.OLS(Data_GMS, sm.add_constant(Data_Precip)).fit()
		Linreg_Coef   = Linreg.params
		Linreg_pvalue = Linreg.pvalues
		Linreg_R2     = Linreg.rsquared

		# Plot
		Plot_Data = {\
			'Data_Precip': Data_Precip, \
			'Data_GMS': Data_GMS, \
			'Linreg_Coef': Linreg_Coef, \
			'Linreg_pvalue': Linreg_pvalue, \
			'Linreg_R2': Linreg_R2, \
		}

		Plot_Config = {\
			'Plot_xlim': Get_xlim(i_Run), \
			'Plot_ylim': Get_ylim(i_Run), \
			'Run': i_Run, \
		}

		Plot_Scatterplot(Plot_Data, Plot_Config)