"""
Plot.EOF_OEMGA.py
==========================
Plot the following charts:
1. Vertical profile for EOF
2. Line chart for PC.
3. Histogram for PC.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mplticker
import sklearn.neighbors
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
	ax.plot(Plot_Data['EOF1'], Prep.Get_Lev_Interp()[1:], color='Black', label='EOF1 (Shallow)')
	ax.plot(Plot_Data['EOF2'], Prep.Get_Lev_Interp()[1:], color='Black', linestyle='dashed', label='EOF2 (Deep)')

	# Plot: reference line
	ax.vlines(0, 0, 1000, colors='Grey', linewidth=1.5, zorder=-2)

	# Configuration
	ax.set_xlabel('Normalized EOF of OMEGA')
	ax.set_xlim([-0.5, 0.5])
	ax.set_xticks(np.arange(-0.5, 0.75, 0.25))
	ax.set_ylabel('Pressure (hPa)')
	ax.set_ylim([975, 200])
	ax.set_yticks([975, 850, 700, 500, 300, 200])
	ax.tick_params(axis='y')

	ax.yaxis.set_major_formatter(mplticker.ScalarFormatter())
	ax.yaxis.get_major_formatter().set_scientific(False)
	ax.yaxis.get_major_formatter().set_useOffset(False)

	# Legend
	ax.legend()
	
	# Plot save
	plt.tight_layout()
	Figure_Path = '../output/Output_Figure/Plot.EOF_OMEGA/'
	Figure_Name = 'Plot.EOF_OMEGA.EOF.{Run}.png'.format(Run=Plot_Config['Run'])
	if (not os.path.exists(Figure_Path)): os.makedirs(Figure_Path)
	plt.savefig(Figure_Path + Figure_Name)
	plt.close('all')

	return

def Plot_Linechart(Plot_Data, Plot_Config):

	# Create matplotlib object
	fig, ax = plt.subplots(figsize=(4, 2), dpi=200)

	# Plot: line charts of PCs
	ax.plot(np.arange(1, 281), Plot_Data['PC1'], color='Black', label='PC1 (Shallow)')
	ax.plot(np.arange(1, 281), Plot_Data['PC2'], color='Black', linestyle='dashed', label='PC2 (Deep)')

	# Plot: reference line
	ax.vlines(0, 0, 1000, colors='Grey', linewidth=1.5)

	# Configuration
	ax.set_xlabel(r'Sample dimension (En$\times$Time) of EOF analysis')
	ax.set_xlim([0, 281])
	ax.set_ylabel('PC')
	ax.set_ylim([0, 0.2])
	ax.set_yticks(np.arange(0, 0.25, 0.05))
	ax.tick_params(axis='y')

	# Legend
	ax.legend()
	
	# Plot save
	plt.tight_layout()
	Figure_Path = '../output/Output_Figure/Plot.EOF_OMEGA/'
	Figure_Name = 'Plot.EOF_OMEGA.PC.{Run}.png'.format(Run=Plot_Config['Run'])
	if (not os.path.exists(Figure_Path)): os.makedirs(Figure_Path)
	plt.savefig(Figure_Path + Figure_Name)
	plt.close('all')

	return

def Plot_Histogram(Plot_Data):

	# Create matplotlib object
	fig = plt.figure(figsize=(5, 5), dpi=300)
	gs = fig.add_gridspec(2, hspace=0)
	ax = gs.subplots(sharex=True, sharey=True)

	# Plot: histogram of PCs
	ax[0].hist(Plot_Data['CTL_PC1_Group1'], bins=150, range=(-0.2, 0.2), density=True, color='Blue', label='Smallest ANO P', alpha=0.15)
	ax[0].hist(Plot_Data['CTL_PC2_Group1'], bins=150, range=(-0.2, 0.2), density=True, color='Blue', alpha=0.15)
	ax[0].plot(Plot_Data['X_Sample'], Plot_Data['CTL_PC1_Group1_KDE'], color='Blue', linewidth=1.5, alpha=1.0, zorder=10)
	ax[0].plot(Plot_Data['X_Sample'], Plot_Data['CTL_PC2_Group1_KDE'], color='Blue', linewidth=1.5, alpha=1.0, zorder=10)

	ax[0].hist(Plot_Data['CTL_PC1_Group2'], bins=150, range=(-0.2, 0.2), density=True, color='Lightcoral', label='Greatest ANO P', alpha=0.15)
	ax[0].hist(Plot_Data['CTL_PC2_Group2'], bins=150, range=(-0.2, 0.2), density=True, color='Lightcoral', alpha=0.15)
	ax[0].plot(Plot_Data['X_Sample'], Plot_Data['CTL_PC1_Group2_KDE'], color='Lightcoral', linewidth=1.5, alpha=1.0, zorder=10)
	ax[0].plot(Plot_Data['X_Sample'], Plot_Data['CTL_PC2_Group2_KDE'], color='Lightcoral', linewidth=1.5, alpha=1.0, zorder=10)

	ax[1].hist(Plot_Data['DEF_PC1_Group1'], bins=150, range=(-0.2, 0.2), density=True, color='Blue', label='Smallest ANO P', alpha=0.15)
	ax[1].hist(Plot_Data['DEF_PC2_Group1'], bins=150, range=(-0.2, 0.2), density=True, color='Blue', alpha=0.15)
	ax[1].plot(Plot_Data['X_Sample'], Plot_Data['DEF_PC1_Group1_KDE'], color='Blue', linewidth=1.5, alpha=1.0, zorder=10)
	ax[1].plot(Plot_Data['X_Sample'], Plot_Data['DEF_PC2_Group1_KDE'], color='Blue', linewidth=1.5, alpha=1.0, zorder=10)

	ax[1].hist(Plot_Data['DEF_PC1_Group2'], bins=150, range=(-0.2, 0.2), density=True, color='Lightcoral', label='Greatest ANO P', alpha=0.15)
	ax[1].hist(Plot_Data['DEF_PC2_Group2'], bins=150, range=(-0.2, 0.2), density=True, color='Lightcoral', alpha=0.15)
	ax[1].plot(Plot_Data['X_Sample'], Plot_Data['DEF_PC1_Group2_KDE'], color='Lightcoral', linewidth=1.5, alpha=1.0, zorder=10)
	ax[1].plot(Plot_Data['X_Sample'], Plot_Data['DEF_PC2_Group2_KDE'], color='Lightcoral', linewidth=1.5, alpha=1.0, zorder=10)

	# Plot: mean line
	ax[0].vlines(np.nanmean(Plot_Data['CTL_PC1_Group1']), 0, 1000, colors='Blue', linestyles='dashed', linewidth=1.0)
	ax[0].vlines(np.nanmean(Plot_Data['CTL_PC2_Group1']), 0, 1000, colors='Blue', linestyles='dashed', linewidth=1.0)
	ax[0].vlines(np.nanmean(Plot_Data['CTL_PC1_Group2']), 0, 1000, colors='Lightcoral', linestyles='dashed', linewidth=1.0)
	ax[0].vlines(np.nanmean(Plot_Data['CTL_PC2_Group2']), 0, 1000, colors='Lightcoral', linestyles='dashed', linewidth=1.0)
	ax[1].vlines(np.nanmean(Plot_Data['DEF_PC1_Group1']), 0, 1000, colors='Blue', linestyles='dashed', linewidth=1.0)
	ax[1].vlines(np.nanmean(Plot_Data['DEF_PC2_Group1']), 0, 1000, colors='Blue', linestyles='dashed', linewidth=1.0)
	ax[1].vlines(np.nanmean(Plot_Data['DEF_PC1_Group2']), 0, 1000, colors='Lightcoral', linestyles='dashed', linewidth=1.0)
	ax[1].vlines(np.nanmean(Plot_Data['DEF_PC2_Group2']), 0, 1000, colors='Lightcoral', linestyles='dashed', linewidth=1.0)

	# Configuration
	ax[1].set_xlabel('PC')
	ax[1].set_xlim([0, 0.20])
	ax[1].set_xticks(np.arange(0, 0.25, 0.05))
	ax[0].set_ylabel('Density')
	ax[1].set_ylabel('Density')
	ax[1].set_ylim([0, 100])
	ax[1].set_yticks(np.arange(0, 100, 20))
	ax[1].tick_params(axis='y')

	# Legend
	ax[0].plot([], [], color='Black', linestyle='dashed', label='Mean of \ndistribution')
	ax[0].legend()

	# Text
	bbox_1 = {'facecolor': 'Blue', 'edgecolor': 'none', 'alpha': 0.1, 'pad': 1}
	bbox_2 = {'facecolor': 'Lightcoral', 'edgecolor': 'none', 'alpha': 0.1, 'pad': 1}
	ax[0].text(np.nanmean(Plot_Data['CTL_PC1_Group1']), 99, '{:.3f}'.format(np.nanmean(Plot_Data['CTL_PC1_Group1'])), color='Blue', fontsize=8, bbox=bbox_1, ha='left', va='top')
	ax[0].text(np.nanmean(Plot_Data['CTL_PC2_Group1']), 99, '{:.3f}'.format(np.nanmean(Plot_Data['CTL_PC2_Group1'])), color='Blue', fontsize=8, bbox=bbox_1, ha='left', va='top')
	ax[0].text(np.nanmean(Plot_Data['CTL_PC1_Group2']), 92, '{:.3f}'.format(np.nanmean(Plot_Data['CTL_PC1_Group2'])), color='Lightcoral', fontsize=8, bbox=bbox_2, ha='left', va='top')
	ax[0].text(np.nanmean(Plot_Data['CTL_PC2_Group2']), 92, '{:.3f}'.format(np.nanmean(Plot_Data['CTL_PC2_Group2'])), color='Lightcoral', fontsize=8, bbox=bbox_2, ha='left', va='top')
	ax[1].text(np.nanmean(Plot_Data['DEF_PC1_Group1']), 99, '{:.3f}'.format(np.nanmean(Plot_Data['DEF_PC1_Group1'])), color='Blue', fontsize=8, bbox=bbox_1, ha='left', va='top')
	ax[1].text(np.nanmean(Plot_Data['DEF_PC2_Group1']), 99, '{:.3f}'.format(np.nanmean(Plot_Data['DEF_PC2_Group1'])), color='Blue', fontsize=8, bbox=bbox_1, ha='left', va='top')
	ax[1].text(np.nanmean(Plot_Data['DEF_PC1_Group2']), 92, '{:.3f}'.format(np.nanmean(Plot_Data['DEF_PC1_Group2'])), color='Lightcoral', fontsize=8, bbox=bbox_2, ha='left', va='top')
	ax[1].text(np.nanmean(Plot_Data['DEF_PC2_Group2']), 92, '{:.3f}'.format(np.nanmean(Plot_Data['DEF_PC2_Group2'])), color='Lightcoral', fontsize=8, bbox=bbox_2, ha='left', va='top')

	# Text: 'CTL' in the left-top corner of the plot (ax[0]) and 'DEF' in the left-top corner of the plot (ax[1])
	ax[0].text(0.01, 0.99, 'CTL', color='Black', fontsize=12, ha='left', va='top', transform=ax[0].transAxes)
	ax[1].text(0.01, 0.99, 'DEF', color='Black', fontsize=12, ha='left', va='top', transform=ax[1].transAxes)
	
	# Plot save
	plt.tight_layout()
	Figure_Path = '../output/Output_Figure/Plot.EOF_OMEGA/'
	Figure_Name = 'Plot.EOF_OMEGA.PC_Histogram.png'
	if (not os.path.exists(Figure_Path)): os.makedirs(Figure_Path)
	plt.savefig(Figure_Path + Figure_Name)
	plt.close('all')

	return

if (__name__ == '__main__'):

	# Get data: EOF and PC
	EOF_CTL = pd.read_csv('../output/Output_Data/EOF_Convection/EOF.CTL.csv')
	PC_CTL = pd.read_csv('../output/Output_Data/EOF_Convection/PC.CTL.csv')
	PC_DEF = pd.read_csv('../output/Output_Data/EOF_Convection/PC.DEF.csv')

	# Plot: vertical profile
	Plot_VProfile(\
		{\
			'EOF1': EOF_CTL['EOF1'], \
			'EOF2': EOF_CTL['EOF2'], \
		}, \
		{\
			'Run': 'CTL', \
		}, \
	)

	# Plot: line chart
	Plot_Linechart(\
		{\
			'PC1': PC_CTL['PC1'], \
			'PC2': PC_CTL['PC2'], \
		}, \
		{\
			'Run': 'CTL', \
		}, \
	)
	Plot_Linechart(\
		{\
			'PC1': PC_DEF['PC1'], \
			'PC2': PC_DEF['PC2'], \
		}, \
		{\
			'Run': 'DEF', \
		}, \
	)

	# Plot: histogram
	Data_PRECT = pd.read_csv('../output/Output_Data/Group/Group.1970-2005.csv')
	X_Sample   = np.arange(0, 0.201, 0.001)
	KDE_BW     = 0.007

	Plot_Data = {}
	Plot_Data['X_Sample']           = X_Sample

	PC1_Group1 = PC_CTL['PC1'][(Data_PRECT['Group']==0)].to_numpy()
	PC1_Group2 = PC_CTL['PC1'][(Data_PRECT['Group']==3)].to_numpy()
	PC2_Group1 = PC_CTL['PC2'][(Data_PRECT['Group']==0)].to_numpy()
	PC2_Group2 = PC_CTL['PC2'][(Data_PRECT['Group']==3)].to_numpy()
	
	Plot_Data['CTL_PC1_Group1']     = PC1_Group1
	Plot_Data['CTL_PC1_Group2']     = PC1_Group2
	Plot_Data['CTL_PC2_Group1']     = PC2_Group1
	Plot_Data['CTL_PC2_Group2']     = PC2_Group2
	Plot_Data['CTL_PC1_Group1_KDE'] = np.exp(sklearn.neighbors.KernelDensity(kernel='gaussian', bandwidth=KDE_BW).fit(PC1_Group1[:, None]).score_samples(X_Sample[:, None]))
	Plot_Data['CTL_PC1_Group2_KDE'] = np.exp(sklearn.neighbors.KernelDensity(kernel='gaussian', bandwidth=KDE_BW).fit(PC1_Group2[:, None]).score_samples(X_Sample[:, None]))
	Plot_Data['CTL_PC2_Group1_KDE'] = np.exp(sklearn.neighbors.KernelDensity(kernel='gaussian', bandwidth=KDE_BW).fit(PC2_Group1[:, None]).score_samples(X_Sample[:, None]))
	Plot_Data['CTL_PC2_Group2_KDE'] = np.exp(sklearn.neighbors.KernelDensity(kernel='gaussian', bandwidth=KDE_BW).fit(PC2_Group2[:, None]).score_samples(X_Sample[:, None]))

	PC1_Group1 = PC_DEF['PC1'][(Data_PRECT['Group']==0)].to_numpy()
	PC1_Group2 = PC_DEF['PC1'][(Data_PRECT['Group']==3)].to_numpy()
	PC2_Group1 = PC_DEF['PC2'][(Data_PRECT['Group']==0)].to_numpy()
	PC2_Group2 = PC_DEF['PC2'][(Data_PRECT['Group']==3)].to_numpy()

	Plot_Data['DEF_PC1_Group1']     = PC1_Group1
	Plot_Data['DEF_PC1_Group2']     = PC1_Group2
	Plot_Data['DEF_PC2_Group1']     = PC2_Group1
	Plot_Data['DEF_PC2_Group2']     = PC2_Group2
	Plot_Data['DEF_PC1_Group1_KDE'] = np.exp(sklearn.neighbors.KernelDensity(kernel='gaussian', bandwidth=KDE_BW).fit(PC1_Group1[:, None]).score_samples(X_Sample[:, None]))
	Plot_Data['DEF_PC1_Group2_KDE'] = np.exp(sklearn.neighbors.KernelDensity(kernel='gaussian', bandwidth=KDE_BW).fit(PC1_Group2[:, None]).score_samples(X_Sample[:, None]))
	Plot_Data['DEF_PC2_Group1_KDE'] = np.exp(sklearn.neighbors.KernelDensity(kernel='gaussian', bandwidth=KDE_BW).fit(PC2_Group1[:, None]).score_samples(X_Sample[:, None]))
	Plot_Data['DEF_PC2_Group2_KDE'] = np.exp(sklearn.neighbors.KernelDensity(kernel='gaussian', bandwidth=KDE_BW).fit(PC2_Group2[:, None]).score_samples(X_Sample[:, None]))

	Plot_Histogram(Plot_Data)