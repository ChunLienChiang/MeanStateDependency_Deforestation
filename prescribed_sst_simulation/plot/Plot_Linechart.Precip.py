"""
Plot_Linechart.Precip.py
========================
Plot line chart of CTL and ANO precipitation
"""

import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import json
import os
import sys
sys.path.append('../')
import preprocessing.Preprocessing_Get_Data as PrepGD
import preprocessing.Preprocessing as Prep

# Get configuration
Config = json.load(open(os.path.join(os.path.dirname(__file__), '../config.json')))

def Plot_Linechart(Plot_Data, Plot_Config):

	# Create matplotlib object
	fig, ax = plt.subplots(figsize=(4, 4), dpi=300)

	# Plot: linechart of precipitation
	ax.plot(\
		np.arange(1, 13), Plot_Data['Precip'], \
		linewidth=2.5, color='Black', \
	)

	# Plot: confidence interval
	ax.fill_between(\
		np.arange(1, 13), *Plot_Data['Precip_CI'], \
		color='Black', ec='none', alpha=0.15, \
	)
	
	# Plot: reference line
	ax.vlines(0, 0, 1000, colors='Grey', linewidth=1.5, zorder=-2)

	# Configuration
	ax.set_xlim([1, 12])
	ax.set_xticks(np.arange(1, 13))
	ax.set_xticklabels(['Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec', 'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun'], rotation=45)
	ax.set_ylabel(r'Precipitation ($W\cdot m^{-2}$)')
	ax.set_ylim(Plot_Config['Plot_ylim'])
	ax.set_yticks(np.linspace(Plot_Config['Plot_ylim'][0], Plot_Config['Plot_ylim'][1], 5))
	ax.tick_params(axis='y')
	
	# Plot save
	plt.tight_layout()
	Figure_Path = '../output/Output_Figure/Plot_Linechart.Precip/'
	Figure_Name = 'Plot_Linechart.Precip.{Run}.png'.format(Run=Plot_Config['Run'])
	if (not os.path.exists(Figure_Path)): os.makedirs(Figure_Path)
	plt.savefig(Figure_Path + Figure_Name)
	plt.close('all')

	return

def Get_ylim(Run):

	if (Run == 'CTL'):

		ylim = [150, 250]

	elif (Run == 'ANO'):

		ylim = [0, 40]

	return ylim

if (__name__ == '__main__'):

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

		# Calculate spatial average
		Data_Precip[i_Run] = Prep.Calc_SpatialAverage(Data_Precip[i_Run], RangeMask='MC_Analysis', Range_Original='MC_Extended_Analysis')

	# Calculate ANO: DEF minus CTL
	Data_Precip['ANO'] = Data_Precip['DEF'] - Data_Precip['CTL']

	# Plot: precipitation
	for i_Run in ['CTL', 'ANO']:

		# Get data
		Plot_Data = {
			'Precip'   : Prep.Calc_CI(Data_Precip[i_Run].reshape((-1, 12)), Test_Axis=0)[0], \
			'Precip_CI': Prep.Calc_CI(Data_Precip[i_Run].reshape((-1, 12)), Test_Axis=0)[1], \
		}

		# Get configuration
		Plot_Config = {
			'Run'      : i_Run, \
			'Plot_ylim': Get_ylim(i_Run), \
		}

		# Plot
		Plot_Linechart(Plot_Data, Plot_Config)