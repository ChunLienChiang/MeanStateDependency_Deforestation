# VProfile_MSE_Budget_Group.py
import numpy as np
import pandas as pd
import scipy.stats as SciSt
import warnings
import Preprocessing as Prep
import Preprocessing_Get_Data as PrepGetData

def Get_XTicks(Var, Run):

	if (Var == 'V_d_MSE'):

		if (Run == 'CTL'):   XTicks = np.arange(-0.015, 0.020, 0.005)
		elif ('ANO' in Run): XTicks = np.arange(-0.005, 0.01, 0.005)
		else: XTicks                = None
	
	elif (Var == 'OMEGA_d_MSE'):

		if (Run == 'CTL'):   XTicks = np.arange(-0.03, 0.04, 0.01)
		elif ('ANO' in Run): XTicks = np.arange(-0.01, 0.015, 0.005)
		else: XTicks                = None
	
	elif (Var == 'V_d_Q'):

		if (Run == 'CTL'):   XTicks = np.arange(-0.05, 0.075, 0.025)
		elif ('ANO' in Run): XTicks = np.arange(-0.015, 0.020, 0.005)
		else: XTicks                = None
	
	elif (Var == 'OMEGA_d_Q'):

		if (Run == 'CTL'):   XTicks = np.arange(-0.03, 0.04, 0.01)
		elif ('ANO' in Run): XTicks = np.arange(-0.01, 0.015, 0.005)
		else: XTicks                = None
	
	elif (Var == 'OMEGA'):

		if (Run == 'CTL'):   XTicks = np.arange(-0.05, 0.075, 0.025)
		elif ('ANO' in Run): XTicks = np.arange(-0.02, 0.03, 0.01)
		else: XTicks                = None
	
	elif (Var == 'd_V'):

		if (Run == 'CTL'):   XTicks = np.arange(-5e-6, 7.5e-6, 2.5e-6)
		elif ('ANO' in Run): XTicks = np.arange(-1e-6, 1.5e-6, 5e-7)
		else: XTicks                = None
	
	elif (Var == 'MSE'):

		if (Run == 'CTL'):   XTicks = np.arange(-400, 600, 200)
		elif ('ANO' in Run): XTicks = np.arange(-400, 600, 200)
		else: XTicks                = None
	
	else: XTicks                    = None

	return XTicks

def Plot_VProfile_Group(Plot_Data, Plot_Config):

	import matplotlib.pyplot as plt
	import matplotlib.ticker as mplticker
	import os

	# Set the color and style of the line plot

	# Create plot objects
	fig, ax = plt.subplots(figsize=(8, 8), dpi=150)

	for ind_Group, i_Group in enumerate(Plot_Data['Data'].keys()):

		i_Line_Color = Plot_Data['Line_Color'][ind_Group]

		# Line plot
		ax.plot(Plot_Data['Data'][i_Group], Prep.Get_Model_Properties()['Lev_Interp'], \
				color=i_Line_Color, marker='.', label=Plot_Config['Legend_Labels'][ind_Group], markersize=20, linewidth=6)

		if not (Plot_Data['Data_CI'][i_Group] is None):

			ax.fill_betweenx(Prep.Get_Model_Properties()['Lev_Interp'], *Plot_Data['Data_CI'][i_Group], fc=i_Line_Color, alpha=0.2, linewidth=6)

	ax.vlines(0, 1000, 0, colors='Grey', linewidth=1.0, zorder=-10)

	if ('YAxis_Signs' in Plot_Config):

		for i_Sig in Plot_Config['YAxis_Signs']:
			
			ax.plot(np.where(Plot_Config['YAxis_Signs'][i_Sig], Plot_Config['Plot_XTicks'][-1]*0.03+Plot_Config['Plot_XTicks'][0]*0.97, np.nan), Prep.Get_Model_Properties()['Lev_Interp'], color='Grey', linestyle='None', alpha=0.9, marker=i_Sig, markersize=10, zorder=15)

	# Plot configuration
	Legend = ax.legend(loc='upper right', prop={'size': 16})
	Legend.set_zorder(50)

	# X axis
	ax.tick_params(axis='x', labelsize=20)
	ax.set_xlabel('{} ({})'.format(Plot_Config['VarLongName'], Plot_Config['VarUnit']), fontsize=20)

	if not (Plot_Config['Plot_XTicks'] is None):

		ax.set_xticks(Plot_Config['Plot_XTicks'])
		ax.set_xlim(Plot_Config['Plot_XTicks'][[0, -1]])

	# Y axis
	ax.set_yscale('log')
	ax.set_yticks(np.concatenate((np.array([975]), np.arange(900, 0, -100))))
	ax.set_yticks(np.concatenate((np.array([975]), np.arange(950, 0, -50))), minor=True)
	ax.tick_params(axis='y', labelsize=20)
	ax.set_ylabel('Pressure (hPa)', fontsize=20)
	ax.set_ylim([975, 400])

	ax.yaxis.set_major_formatter(mplticker.ScalarFormatter())
	ax.yaxis.get_major_formatter().set_scientific(False)
	ax.yaxis.get_major_formatter().set_useOffset(False)
	ax.yaxis.set_minor_formatter(mplticker.NullFormatter())

	ax.minorticks_on()
	ax.grid(which='major', color='grey', linewidth=0.6, alpha=0.6, zorder=15)
	ax.grid(which='minor', color='grey', linewidth=0.3, alpha=0.4, zorder=15)
	
	# Save figure
	Figure_Path = 'Output_Figures/VProfile_MSE_Budget_Group/'
	Figure_Name = 'VProfile_MSE_Budget_Group_{}_{}_{}.png'.format(Plot_Config['Plot_Month'], Plot_Config['Var'], Plot_Config['Run'])

	if not (os.path.exists(Figure_Path)):

		os.mkdir(Figure_Path)

	plt.tight_layout()
	plt.savefig(Figure_Path + Figure_Name, bbox_inches='tight')
	plt.close('all')

	return

def Get_Data(Run, Var, Events_List, MiPert_List=np.arange(5)):

	Data = []

	for i_Event in Events_List:

		for i_MiPert in MiPert_List:
			
			RawData = PrepGetData.Get_MSE_Budget(Run, i_Event, Var, MiPert=i_MiPert)
			Data.append(RawData)
	
	return Data

def Calc_SigLevel(Mean_1, Std_1, n_1, Mean_2, Std_2, n_2):

	Diff_Pvalue = SciSt.ttest_ind_from_stats(Mean_1, Std_1, n_1, Mean_2, Std_2, n_2, equal_var=False, alternative='greater').pvalue
	SignificanceLevel_095 = np.where(Diff_Pvalue<0.05, True, False)
	SignificanceLevel_090 = np.where(Diff_Pvalue<0.10, True, False)

	return {'<': np.where(SignificanceLevel_095, True, False), 'o': np.where(SignificanceLevel_090&(~SignificanceLevel_095), True, False)}

if (__name__ == '__main__'):

	warnings.simplefilter('ignore', category=RuntimeWarning)
	
	Model_Properties = Prep.Get_Model_Properties()

	# Get data
	Var_List = [\
		'V_d_MSE', \
		'OMEGA_d_MSE', \
		'MSE', \
		'OMEGA', \
		'd_V', \
	]
	VarLongName_List = [\
		'Horizontal advection of MSE', \
		'Vertical advection of MSE', \
		'MSE', \
		'Vertical wind velocity', \
		'Convergence of horizontal wind', \
	]
	VarUnit_List = [\
		'$W\cdot kg^{-1}$', \
		'$W\cdot kg^{-1}$', \
		'$J\cdot kg^{-1}$', \
		'$Pa\cdot s^{-1}$', \
		'$s^{-1}$', \
	]

	
	n_Group     = 4
	Plot_Group  = [0, 3]
	'''
	n_Group     = 4
	Plot_Group  = [-1, 1]
	'''
	MiPert_List  = np.array([0])
	n_MiPert     = np.size(MiPert_List)

	# Get the list of events
	Events_List = PrepGetData.Get_Events_List()

	for i_Month in ['All', 'DJF', 'JA']:

		# Month mask
		Mask_Month        = Prep.Get_Mask_Month(i_Month, ENSO_Shifting=True)

		# ANO P group
		EventsGroup       = pd.read_csv('Output_Data/File_EventsGroup/File_EventsGroup_{}Gs.csv'.format(n_Group))
		EventsGroup_List  = EventsGroup.loc[EventsGroup['Month']==i_Month, 'Group_ANOP_ANOSF'].to_numpy()
		#EventsGroup_List  = EventsGroup.loc[EventsGroup['Month']==i_Month, 'Group_Nino34_ANOSF'].to_numpy()

		# ANO P list
		ANOP_List         = pd.read_csv('Output_Data/File_PRECT/File_PRECT.csv')
		ANOP_List         = ANOP_List.loc[(ANOP_List['Compset']=='piControl')&(ANOP_List['Month']==i_Month), 'ANOSF'].to_numpy()

		Data = {}
	
		# Get data
		for i_Var, i_VarLongName, i_VarUnit in zip(Var_List, VarLongName_List, VarUnit_List):
			
			#################################
			# Plot for transient simulation #

			Data[i_Var] = {}

			for i_Run in ['CTL', 'DEFSF']:

				Data[i_Var][i_Run] = Get_Data(i_Run, i_Var, Events_List, MiPert_List=MiPert_List)
				Data[i_Var][i_Run] = [np.nanmean(np.where(Mask_Month[:, None], i, np.nan), axis=0) for i in Data[i_Var][i_Run]]
				Data[i_Var][i_Run] = np.concatenate(tuple([i[None, :] for i in Data[i_Var][i_Run]]), axis=0)
				
			for i_Run in ['CTL', 'ANOSF']:

				Plot_Data_Mean = {}
				Plot_Data_CI   = {}
				Plot_Data_Std  = {}
				Plot_Data_n    = {}
				Plot_ANOP      = {}

				for i_Group in Plot_Group:
					
					if (i_Run == 'CTL'):
						
						if (i_Var== 'MSE'):

							Data_Masked = np.where(np.isin(EventsGroup_List[:, None], i_Group), Data[i_Var][i_Run], np.nan) - np.nanmean(Data[i_Var][i_Run], axis=0)
						
						else:
							
							Data_Masked = np.where(np.isin(EventsGroup_List[:, None], i_Group), Data[i_Var][i_Run], np.nan)
						
						Data_Mean, Data_CI, Data_Std, Data_n = Prep.Calc_CI(Data_Masked, Bool_return_Std=True)
					
					elif (i_Run == 'ANOSF'):

						Data_Masked_1 = np.where(np.isin(EventsGroup_List[:, None], i_Group), Data[i_Var]['DEFSF'], np.nan)
						Data_Masked_2 = np.where(np.isin(EventsGroup_List[:, None], i_Group), Data[i_Var]['CTL'], np.nan)
						Data_Mean, Data_CI, Data_Std, Data_n = Prep.Calc_CI([Data_Masked_1, Data_Masked_2], Bool_return_Std=True)
					
					else:
						
						Data_Mean = np.where(np.isin(EventsGroup_List[:, None], i_Group), Data[i_Var][i_Run], np.nan)
						Data_CI   = None
					
					if (isinstance(i_Group, list)):

						Plot_Data_Mean['Group_{}/{}'.format(*i_Group)] = Data_Mean
						Plot_Data_CI['Group_{}/{}'.format(*i_Group)]   = Data_CI
						Plot_Data_Std['Group_{}/{}'.format(*i_Group)]  = Data_Std
						Plot_Data_n['Group_{}/{}'.format(*i_Group)]    = Data_n
						Plot_ANOP['Group_{}/{}'.format(*i_Group)]      = np.nanmean(np.where(np.isin(EventsGroup_List, i_Group), ANOP_List, np.nan))

					else:

						Plot_Data_Mean['Group_{}'.format(int(i_Group))] = Data_Mean
						Plot_Data_CI['Group_{}'.format(int(i_Group))]   = Data_CI
						Plot_Data_Std['Group_{}'.format(int(i_Group))]  = Data_Std
						Plot_Data_n['Group_{}'.format(int(i_Group))]    = Data_n
						Plot_ANOP['Group_{}'.format(int(i_Group))]      = np.nanmean(np.where(np.isin(EventsGroup_List, i_Group), ANOP_List, np.nan))

				# Calculate significance level
				YAxis_Signs = Calc_SigLevel(\
					Plot_Data_Mean['Group_{}'.format(Plot_Group[0])], \
					Plot_Data_Std['Group_{}'.format(Plot_Group[0])], \
					Plot_Data_n['Group_{}'.format(Plot_Group[0])], \
					Plot_Data_Mean['Group_{}'.format(Plot_Group[-1])], \
					Plot_Data_Std['Group_{}'.format(Plot_Group[-1])], \
					Plot_Data_n['Group_{}'.format(Plot_Group[-1])]\
				)

				# Plot
				Plot_Data = {\
					'Data': Plot_Data_Mean, \
					'Data_CI': Plot_Data_CI, \
					'Line_Color': ['C0', 'C3'], \
				}

				Plot_Config = {\
					'Var'          : i_Var, \
					'VarLongName'  : i_VarLongName, \
					'VarUnit'      : i_VarUnit, \
					'Run'          : i_Run, \
					'Plot_Month'   : i_Month, \
					'Plot_XTicks'  : Get_XTicks(i_Var, i_Run), \
					'YAxis_Signs'  : YAxis_Signs, \
					'Legend_Labels': [\
						str(np.round(Plot_ANOP['Group_0'], 1)) + r' $W\cdot m^{-2}$' + ' (pct.0~25), n=' + str(int(np.nanmin(Plot_Data_n['Group_0']))), str(np.round(Plot_ANOP['Group_3'], 1)) + r' $W\cdot m^{-2}$' + ' (pct.75~100), n=' + str(int(np.nanmin(Plot_Data_n['Group_3'])))] \
				}
				
				Plot_VProfile_Group(Plot_Data, Plot_Config)
			