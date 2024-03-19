"""
tool_Group.py
==========================
Calculate the grouping for analysis.
"""

import numpy as np
import pandas as pd
import json
import os
import sys
sys.path.append('../../')
import preprocessing.Preprocessing_Get_Data as PrepGD
import preprocessing.Preprocessing as Prep

# Get configuration
Config = json.load(open(os.path.join(os.path.dirname(__file__), '../../config.json')))

if (__name__ == '__main__'):

	# Set output path
	File_Path = '../../output/Output_Data/Group/'
	File_Name = 'Group.{Year_Start}-{Year_End}.csv'.format(\
		Year_Start=Config['Data_TimeRange'][0], \
		Year_End=Config['Data_TimeRange'][1], \
	)

	if not (os.path.exists(File_Path)): os.makedirs(File_Path)

	# Create empty list
	PRECT_CTL = []
	PRECT_DEF = []
	PRECT_ANO = []

	for i_En in Config['Data_En']:
		
		# Get data: PRECT
		PRECT_CTL_Temp = PrepGD.Get_Simulation_Extract('CTL', i_En, 'PRECT', 'MC_Extended_Analysis')
		PRECT_CTL_Temp = Prep.Calc_SpatialAverage(PRECT_CTL_Temp, LandMask='Land', RangeMask='MC_Analysis', Range_Original='MC_Extended_Analysis')[6:-6, ...]
		PRECT_CTL_Temp = Prep.Calc_AnnualMean(PRECT_CTL_Temp)
		
		PRECT_DEF_Temp = PrepGD.Get_Simulation_Extract('DEF', i_En, 'PRECT', 'MC_Extended_Analysis')
		PRECT_DEF_Temp = Prep.Calc_SpatialAverage(PRECT_DEF_Temp, LandMask='Land', RangeMask='MC_Analysis', Range_Original='MC_Extended_Analysis')[6:-6, ...]
		PRECT_DEF_Temp = Prep.Calc_AnnualMean(PRECT_DEF_Temp)

		PRECT_CTL.append(PRECT_CTL_Temp)
		PRECT_DEF.append(PRECT_DEF_Temp)
		PRECT_ANO.append(PRECT_DEF_Temp - PRECT_CTL_Temp)

	PRECT_CTL = np.concatenate(tuple([i[None, :] for i in PRECT_CTL]), axis=0).ravel()
	PRECT_DEF = np.concatenate(tuple([i[None, :] for i in PRECT_DEF]), axis=0).ravel()
	PRECT_ANO = np.concatenate(tuple([i[None, :] for i in PRECT_ANO]), axis=0).ravel()

	# Calculate grouping
	Group   = np.full_like(PRECT_ANO, np.nan)
	n_Group = 4

	for i_Group in np.arange(n_Group):

		Group[(PRECT_ANO>=np.nanpercentile(PRECT_ANO, (100/n_Group)*i_Group))&(PRECT_ANO<=np.nanpercentile(PRECT_ANO, (100/n_Group)*(i_Group+1)))] = i_Group

	# Output to csv file
	Grouping = pd.DataFrame({\
		'En'   : np.repeat(np.array(Config['Data_En']), PRECT_ANO.size/len(Config['Data_En'])), \
		'Year' : np.tile(np.arange(*Config['Data_TimeRange']), len(Config['Data_En'])), \
		'Group': Group.astype(int), \
		'PRECT_CTL': PRECT_CTL, \
		'PRECT_DEF': PRECT_DEF, \
		'PRECT_ANO': PRECT_ANO, \
	})
	Grouping.to_csv(File_Path + File_Name, index=False)