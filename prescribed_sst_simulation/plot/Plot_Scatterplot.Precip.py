"""
Plot_Scatterplot.Precip.py
==========================
Plot scatterplot between CTL and ANOprecipitation.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import statsmodels.api as sm
import os
import json
import sys
sys.path.append('../')
import preprocessing.Preprocessing_Get_Data as PrepGD
import preprocessing.Preprocessing as Prep

# Get configuration
Config = json.load(open(os.path.join(os.path.dirname(__file__), '../config.json')))

def Plot_Scatterplot(Plot_Data, Plot_Config):

	# Create matplotlib object
	fig, ax = plt.subplots(figsize=(4, 4), dpi=300)

	# Plot: scatterplot
	ax.scatter(\
		Plot_Data['CTL'], \
		Plot_Data['ANO'], \
		s=10, \
		color='Blue', \
		alpha=1.0, \
	)

	# Plot: regression line
	ax.plot(\
		np.array(Plot_Config['Plot_xticks'][[0, -1]]), \
		Plot_Data['Linreg_Coef'][0] + Plot_Data['Linreg_Coef'][1]*np.array(Plot_Config['Plot_xticks'][[0, -1]]), \
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
	ax.set_xlabel(r'CTL Precipitation ($W\cdot m^{-2}$)')
	ax.set_xlim(Plot_Config['Plot_xticks'][[0, -1]])
	ax.set_xticks(Plot_Config['Plot_xticks'])
	ax.set_ylabel(r'ANO Precipitation ($W\cdot m^{-2}$)')
	ax.set_ylim(Plot_Config['Plot_yticks'][[0, -1]])
	ax.set_yticks(Plot_Config['Plot_yticks'])

	# Plot save
	plt.tight_layout()
	Figure_Path = '../output/Output_Figure/Plot_Scatterplot.Precip/'
	if (Plot_Config['Season'] is None): Figure_Name = 'Plot_Scatterplot.Precip.png'
	else: Figure_Name = 'Plot_Scatterplot.Precip.{Season}.png'.format(Season=Plot_Config['Season'])
	if (not os.path.exists(Figure_Path)): os.makedirs(Figure_Path)
	plt.savefig(Figure_Path + Figure_Name)
	plt.close('all')

	return

def Get_xticks(Season):

	if (Season is None):
		
		xticks = np.arange(125, 325, 50)
	
	else:
		
		xticks = np.arange(0, 400, 100)

	return xticks

def Get_yticks(Season):

	if (Season is None):
		
		yticks = np.arange(-25, 75, 25)
	
	else:
		
		yticks = np.arange(-50, 150, 50)

	return yticks

if (__name__ == '__main__'):

	# Plot for different seasons
	for i_Season in [None, 'DJF', 'JA']:
		
		# Read data: Precipitation
		Data_Precip = {}
		
		for i_Run in ['CTL', 'DEF']:

			Data_Precip[i_Run] = []

			for i_En in Config['Data_En']:

				Data_Precip[i_Run].append(PrepGD.Get_Simulation_Extract(\
					i_Run, \
					i_En, \
					'PRECT', \
					'MC_Extended_Analysis', \
				))
			
			# Concatenate and remove the first 6 months and the last 6 months
			Data_Precip[i_Run] = np.concatenate(tuple([i[6:-6, ...] for i in Data_Precip[i_Run]]), axis=0)

			# Calculate annual mean
			Data_Precip[i_Run] = Prep.Calc_AnnualMean(Data_Precip[i_Run], Season=i_Season)

			# Calculate spatial average
			Data_Precip[i_Run] = Prep.Calc_SpatialAverage(Data_Precip[i_Run], RangeMask='MC_Analysis', Range_Original='MC_Extended_Analysis')
		
		# Calculate ANO: DEF minus CTL
		Data_Precip['ANO'] = Data_Precip['DEF'] - Data_Precip['CTL']

		# Calculate regression line
		Linreg        = sm.OLS(Data_Precip['ANO'], sm.add_constant(Data_Precip['CTL'])).fit()
		Linreg_Coef   = Linreg.params
		Linreg_pvalue = Linreg.pvalues
		Linreg_R2     = Linreg.rsquared
		
		# Plot: scatterplot
		Plot_Data = {\
			'CTL': Data_Precip['CTL'], \
			'ANO': Data_Precip['ANO'], \
			'Linreg_Coef': Linreg_Coef, \
			'Linreg_pvalue': Linreg_pvalue, \
			'Linreg_R2': Linreg_R2, \
		}

		Plot_Config = {\
			'Plot_xticks': Get_xticks(i_Season), \
			'Plot_yticks': Get_yticks(i_Season), \
			'Season': i_Season, \
		}

		Plot_Scatterplot(\
			Plot_Data, \
			Plot_Config, \
		)