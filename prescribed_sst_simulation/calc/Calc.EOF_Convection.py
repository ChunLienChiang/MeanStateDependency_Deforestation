"""
Calc.EOF_Convection.py
==========================
Calculate the EOF of CTL omega and extract modes of shallow convection and deep convection.
Project DEF omega onto CTL EOF and get pseudo PC.
"""

import numpy as np
import pandas as pd
import eofs.standard
import statsmodels.api as sm
import os
import json
import sys
sys.path.append('../')
import preprocessing.Preprocessing_Get_Data as PrepGD
import preprocessing.Preprocessing as Prep

# Get configuration
Config = json.load(open(os.path.join(os.path.dirname(__file__), '../config.json')))

def Calc_EOF(Data):

	EOF_Solver = eofs.standard.Eof(Data, center=False)
	EOF        = EOF_Solver.eofs(neofs=2, eofscaling=0).filled(np.nan)
	PC         = EOF_Solver.pcs(npcs=2, pcscaling=0)

	return EOF, PC

def Calc_EOF_Transform(e1, e2, pc1, pc2):

	"""
	Linear transformation: Back and Bretherton (2009)
	"""

	D1         = e1[0]
	D2         = e2[0]
	EOF2_Trans = (-D2/(D1-D2)) * e1 + (D1/(D1-D2)) * e2
	V1         = np.dot(e1, EOF2_Trans)
	V2         = np.dot(e2, EOF2_Trans)
	EOF1_Trans = ((-V2*(D1+D2))/(V1*D2-V2*D1)) * e1 + ((V1*(D1+D2))/(V1*D2-V2*D1)) * e2
	PC1_Trans  = ((-V2*(D1+D2))/(V1*D2-V2*D1)) * pc1 + ((V1*(D1+D2))/(V1*D2-V2*D1)) * pc2
	PC2_Trans  = (-D2/(D1-D2)) * pc1 + (D1/(D1-D2)) * pc2
	
	return EOF1_Trans, EOF2_Trans, PC1_Trans, PC2_Trans

def Print_OLS_Summary(X, Y):
	
	Model_Results = sm.OLS(Y, sm.add_constant(X)).fit()
	
	File_Path = '../output/Output_Data/EOF_Convection/'
	print(Model_Results.summary(), file=open(File_Path+'Linreg_OLS_Summary.txt', 'w'))

	return Model_Results

def Output_csv(EOF, PC, Run):

	File_Path = '../output/Output_Data/EOF_Convection/'
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

	# Get data
	OMEGA = {\
		'CTL': [], \
		'DEF': [], \
	}

	for i_Run in ['CTL', 'DEF']:

		for i_En in Config['Data_En']:

			# Get data: Omega
			OMEGA_Temp = PrepGD.Get_Simulation_Extract(i_Run, i_En, 'OMEGA', 'MC_Extended_Analysis')
			OMEGA_Temp = Prep.Calc_SpatialAverage(OMEGA_Temp, LandMask='Land', RangeMask='MC_Analysis', Range_Original='MC_Extended_Analysis')
			OMEGA_Temp = Prep.Calc_AnnualMean(OMEGA_Temp[6:-6, ...], 0)
			OMEGA[i_Run].append(OMEGA_Temp[..., 1:])
	
		OMEGA[i_Run] = np.concatenate(OMEGA[i_Run], axis=0)

	# Calculate EOF
	EOF_CTL, PC_CTL = Calc_EOF(OMEGA['CTL'])
	
	# Linear transformation
	EOF1_Trans, EOF2_Trans, PC1_Trans, PC2_Trans = Calc_EOF_Transform(EOF_CTL[0, ...], EOF_CTL[1, ...], PC_CTL[..., 0], PC_CTL[..., 1])
	EOF_CTL = np.concatenate(tuple([EOF1_Trans[None, :], EOF2_Trans[None, :]]), axis=0)
	
	# Change sign of EOF and PC if the EOF (vertical profile of omega) between 800 hPa to 975 hPa is positive
	for i_PC in np.arange(2):
		
		if (np.nanmean(EOF_CTL[i_PC, :8])>=0):

			PC_CTL[..., i_PC]  = PC_CTL[..., i_PC] * -1
			EOF_CTL[i_PC, ...] = EOF_CTL[i_PC, ...] * -1
	
	# EOF projection for DEF omega to get pseudo-PC
	PC_CTL = np.dot(OMEGA['CTL'], EOF_CTL.T)
	PC_DEF = np.dot(OMEGA['DEF'], EOF_CTL.T)

	# Save EOF results
	Output_csv(EOF_CTL, PC_CTL, 'CTL')
	Output_csv(None, PC_DEF, 'DEF')

	# Calculate OLS model
	Data_PRECT = pd.read_csv('../output/Output_Data/Group/Group.1970-2005.csv')
	
	PC_CTL_Group0 = pd.DataFrame({'PC1': PC_CTL[..., 0], 'PC2': PC_CTL[..., 1]})[(Data_PRECT['Group']==0)].mean(axis=0).round(3)
	PC_CTL_Group3 = pd.DataFrame({'PC1': PC_CTL[..., 0], 'PC2': PC_CTL[..., 1]})[(Data_PRECT['Group']==3)].mean(axis=0).round(3)
	PC_DEF_Group0 = pd.DataFrame({'PC1': PC_DEF[..., 0], 'PC2': PC_DEF[..., 1]})[(Data_PRECT['Group']==0)].mean(axis=0).round(3)
	PC_DEF_Group3 = pd.DataFrame({'PC1': PC_DEF[..., 0], 'PC2': PC_DEF[..., 1]})[(Data_PRECT['Group']==3)].mean(axis=0).round(3)
	print(PC_CTL_Group0)
	print(PC_CTL_Group3)
	print(PC_DEF_Group0)
	print(PC_DEF_Group3)
	print(PC_DEF_Group0 - PC_CTL_Group0)
	print(PC_DEF_Group3 - PC_CTL_Group3)

	Model_Results = Print_OLS_Summary(\
		pd.DataFrame({'PC1': PC_CTL[..., 0], 'PC2': PC_CTL[..., 1]}), \
		Data_PRECT['PRECT_CTL'], \
	)