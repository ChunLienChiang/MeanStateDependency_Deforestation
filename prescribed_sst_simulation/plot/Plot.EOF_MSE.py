"""
Plot.EOF_MSE.py
==========================
Plot the following charts:
1. Vertical profile for EOF
2. Line chart for PC.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mplticker
import json
import os
import sys
sys.path.append('../')
import preprocessing.Preprocessing as Prep

# Get configuration
Config = json.load(open(os.path.join(os.path.dirname(__file__), '../config.json')))

def Plot_VProfile(Plot_Data, Plot_Config):

	# Create matplotlib object
	fig, ax = plt.subplots(figsize=(4, 4), dpi=200)

	# Plot: vertical profiles of EOFs
	ax.plot(Plot_Data, Prep.Get_Lev_Interp(), color='Blue')

	# Plot: reference line
	ax.vlines(0, 0, 1000, colors='Grey', linewidth=1.5, zorder=-2)

	# Configuration
	ax.set_xlabel('EOF')
	ax.set_xlim(Plot_Config['Plot_xlim'])
	ax.set_xticks([Plot_Config['Plot_xlim'][0], 0, Plot_Config['Plot_xlim'][1]])
	ax.set_ylabel('Pressure (hPa)')
	ax.set_ylim([975, 200])
	ax.set_yticks([975, 850, 700, 500, 300, 200])
	ax.tick_params(axis='y')

	ax.yaxis.set_major_formatter(mplticker.ScalarFormatter())
	ax.yaxis.get_major_formatter().set_scientific(False)
	ax.yaxis.get_major_formatter().set_useOffset(False)

	# Plot save
	plt.tight_layout()
	Figure_Path = '../output/Output_Figure/Plot.EOF_MSE/'
	Figure_Name = 'Plot.EOF_MSE.EOF.{Run}.PC{PC}.png'.format(Run=Plot_Config['Run'], PC=Plot_Config['PC'])
	if (not os.path.exists(Figure_Path)): os.makedirs(Figure_Path)
	plt.savefig(Figure_Path + Figure_Name)
	plt.close('all')

	return

def Plot_Linechart(Plot_Data, Plot_Config):

	# Create matplotlib object
	fig, ax = plt.subplots(figsize=(4, 4), dpi=200)

	# Plot: line charts of PCs
	ax.plot(np.arange(*Config['Data_TimeRange']), Plot_Data, color='Blue')

	# Plot: reference line
	ax.vlines(0, 0, 1000, colors='Grey', linewidth=1.5)

	# Configuration
	ax.set_xlabel('Year')
	ax.set_xlim(Config['Data_TimeRange'])
	ax.set_ylabel('PC')
	ax.set_ylim([-2.5, 2.5])
	ax.set_yticks(np.arange(-2.5, 3.0, 0.5))
	ax.tick_params(axis='y')

	# Plot save
	plt.tight_layout()
	Figure_Path = '../output/Output_Figure/Plot.EOF_MSE/'
	Figure_Name = 'Plot.EOF_MSE.PC.{Run}.PC{PC}.png'.format(Run=Plot_Config['Run'], PC=Plot_Config['PC'])
	if (not os.path.exists(Figure_Path)): os.makedirs(Figure_Path)
	plt.savefig(Figure_Path + Figure_Name)
	plt.close('all')

	return

def Get_xlim(Run):

	if (Run == 'CTL'):

		return [-1000, 1000]

	elif (Run == 'ANO'):

		return [-200, 200]

	else:

		raise ValueError('The value of Run is not valid.')

	return

if (__name__ == '__main__'):

	for i_Run in ['CTL', 'ANO']:

		# Get data: EOF and PC
		EOF = pd.read_csv('../output/Output_Data/EOF_MSE/EOF.{Run}.csv'.format(Run=i_Run))
		PC = pd.read_csv('../output/Output_Data/EOF_MSE/PC.{Run}.csv'.format(Run=i_Run))

		for i_PC in np.arange(2):

			# Plot: vertical profile
			Plot_VProfile(\
				EOF['EOF{}'.format(i_PC+1)], \
				{\
					'Plot_xlim': Get_xlim(i_Run), \
					'PC' : i_PC+1, \
					'Run': i_Run, \
				}, \
			)

			# Plot: line chart
			Plot_Linechart(\
				 PC['PC{}'.format(i_PC+1)], \
				{\
					'PC' : i_PC+1, \
					'Run': i_Run, \
				}, \
			)