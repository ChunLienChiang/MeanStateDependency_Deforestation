"""
Plot_Barchart.GMS.py
==========================
Barchart for mean GMS in two groups (the greater/smaller ANO precipitation)
"""

import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
import json
import os
import sys
sys.path.append('../')
import preprocessing.Preprocessing as Prep

# Get configuration
Config = json.load(open(os.path.join(os.path.dirname(__file__), '../config.json')))

def Plot(Plot_Data):

	X_Center = np.array([1, 3, 5])

	# Create matplotlib object
	fig, ax = plt.subplots(figsize=(5, 5), dpi=200)

	# Plot: barchart for blue group
	Bar_1 = ax.bar(\
		X_Center - 0.3, \
		[Plot_Data[Var]['Weakest'][0] for Var in list(Plot_Data.keys())], \
		width=0.45, \
		color=['Blue', 'Blue', 'Blue'], \
		yerr=[(Plot_Data[Var]['Weakest'][1][1]-Plot_Data[Var]['Weakest'][1][0])/2 for Var in list(Plot_Data.keys())], \
		alpha=0.7, \
	)
	
	# Plot: barchart for red group
	Bar_2 = ax.bar(\
		X_Center + 0.3, \
		[Plot_Data[Var]['Greatest'][0] for Var in list(Plot_Data.keys())], \
		width=0.45, \
		color=['Red', 'Red', 'Red'], \
		yerr=[(Plot_Data[Var]['Greatest'][1][1]-Plot_Data[Var]['Greatest'][1][0])/2 for Var in list(Plot_Data.keys())], \
		alpha=0.7, \
	)
	
	# Plot: reference line
	ax.hlines(0, 0, 10, colors='Grey', linewidth=1.5)

	# Legend
	ax.legend((Bar_1[0], Bar_2[0]), ('Smallest ANO P', 'Greatest ANO P'))
	
	# Configuration
	ax.set_xlim([0, 6])
	ax.set_xticks(X_Center)
	ax.set_xticklabels([r'$M^{\prime}$', r'$M^{\prime}_{h_m}$', r'$M^{\prime}_{\Omega}$'])
	ax.set_ylim([-0.01, 0.01])
	ax.set_yticks(np.arange(-0.01, 0.015, 0.005))
	ax.set_ylabel(r'$J\cdot kg^{-1}\cdot Pa^{-1}$')

	# Plot save
	plt.tight_layout()
	Figure_Path = '../output/Output_Figure/Plot_Barchart.GMS/'
	Figure_Name = 'Plot_Barchart.GMS.png'
	if (not os.path.exists(Figure_Path)): os.makedirs(Figure_Path)
	plt.savefig(Figure_Path + Figure_Name)
	plt.close('all')

	return

if (__name__ == '__main__'):
	
	# Get grouping
	Group = pd.read_csv('../output/Output_Data/Group/Group.{Year_Start}-{Year_End}.csv'.format(\
		Year_Start=Config['Data_TimeRange'][0], \
		Year_End=Config['Data_TimeRange'][1], \
	))['Group'].to_numpy()
	
	# Get GMS
	GMS = {\
		'ANO': [], \
		'ANO_ThDyn': [], \
		'ANO_Dyn': [], \
	}

	for i_En in Config['Data_En']:
		
		RawData = xr.open_dataset('../output/Output_Data/GMS/GMS.1970-2005.En{En}.MC_Extended.nc'.format(En=i_En), engine='netcdf4')
		GMS['ANO'].append(RawData['GMS_ANO'].to_numpy())
		GMS['ANO_ThDyn'].append(RawData['GMS_ANO_ThDyn'].to_numpy())
		GMS['ANO_Dyn'].append(RawData['GMS_ANO_Dyn'].to_numpy())
	
	GMS = {key: np.concatenate(GMS[key]) for key in GMS}
	
	# Plot
	Plot_Data = {\
		'GMS_{Run}'.format(Run=Run): {\
			'Weakest': Prep.Calc_CI(np.where(Group==0, GMS[Run], np.nan), Test_Axis=0), \
			'Greatest': Prep.Calc_CI(np.where(Group==3, GMS[Run], np.nan), Test_Axis=0), \
		} for Run in list(GMS.keys())
	}

	Plot(Plot_Data)