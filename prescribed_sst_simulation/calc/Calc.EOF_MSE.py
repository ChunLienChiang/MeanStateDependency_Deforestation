"""
Calc.EOF_MSE.py
==========================
Calculate the EOF of CTL MSE profile.
"""

import numpy as np
import pandas as pd
import xarray as xr
import eofs.standard
import os
import json

# Get configuration
Config = json.load(open(os.path.join(os.path.dirname(__file__), '../config.json')))

def Calc_EOF(Data):

	EOF_Solver = eofs.standard.Eof(Data, center=False)
	EOF        = EOF_Solver.eofs(neofs=2, eofscaling=2)
	PC         = EOF_Solver.pcs(npcs=2, pcscaling=1)

	return EOF, PC

def Output_csv(EOF, PC, Run):

	File_Path = '../output/Output_Data/EOF_MSE/'
	if (not os.path.exists(File_Path)): os.makedirs(File_Path)
	
	if (not EOF is None):
		
		pd.DataFrame({\
			'EOF{}'.format(i+1): EOF[i, ...] for i in np.arange(2) \
		}).to_csv(File_Path + 'EOF.{Run}.csv'.format(Run=Run), index=False)
	
	if (not PC is None):

		pd.DataFrame({\
			'PC{}'.format(i+1): PC[..., i] for i in np.arange(2) \
		}).to_csv(File_Path + 'PC.{Run}.csv'.format(Run=Run), index=False)

	return

if (__name__ == '__main__'):

	for i_Run in ['CTL', 'ANO']:

		# Read data: MSE
		Data_MSE = []

		for i_En in Config['Data_En']:

			Data_File_Path = '../output/Output_Data/MSE_Budget/'
			Data_File_Name = 'MSE_Budget.1970-2005.En{En}.MC_Extended.nc'.format(En=i_En)
			if not (os.path.exists(Data_File_Path + Data_File_Name)): raise ValueError('The file {} does not exist.'.format(Data_File_Name))

			Data_MSE.append(xr.open_dataset(Data_File_Path + Data_File_Name)['MSE_{Run}'.format(Run=i_Run)].to_numpy()[None, ...])

		Data_MSE = np.nanmean(np.concatenate(tuple(Data_MSE), axis=0), axis=0)
		
		# Calculate anomalies of MSE vertical profile
		Data_MSE = Data_MSE - np.nanmean(Data_MSE, axis=0)

		# Calculate EOF
		EOF, PC = Calc_EOF(Data_MSE)

		# Save data
		Output_csv(EOF, PC, i_Run)