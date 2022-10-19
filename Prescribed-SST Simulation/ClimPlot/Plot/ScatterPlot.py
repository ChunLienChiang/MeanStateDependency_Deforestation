# ScatterPlot.py

# Import modules
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

def Figure_Initialization(Plot_Config):

	# Create new figure
	FigureObject, AxesObject = plt.subplots(figsize=(8, 8), dpi=300)

	return FigureObject, AxesObject, Plot_Config

def Calc_LinearPredict(Y, X):

	import statsmodels.api as sm

	Reg_X  = sm.add_constant(X)
	Model  = sm.OLS(Y.astype(float), Reg_X.astype(float))
	Result = Model.fit()

	return Result.predict(sm.add_constant([np.min(X), np.max(X)])), Result.params, Result.pvalues, Result.rsquared

def Draw_Plot(ax, Data, Data_Config, Plot_Config, **kwargs):

	# Plot
	if ('Plot_Cmap' in Plot_Config):

		import matplotlib.colors

		cmap = Plot_Config['Plot_Cmap']
		norm = matplotlib.colors.Normalize(vmin=np.min(Data['Color']), vmax=np.max(Data['Color']))

		New_PlotObject = ax.scatter(Data['X'], Data['Y'], c=Data['Color'], cmap=cmap, norm=norm, s=40, alpha=0.3, zorder=10)
		
		if ('X_Primary' in Data):

			New_PlotObject = ax.scatter(Data['X_Primary'], Data['Y_Primary'], c=Data['Color_Primary'], cmap=cmap, norm=norm, \
										s=100, edgecolors='Black', linewidths=3, alpha=0.9, zorder=15)
		
	else:

		New_PlotObject = ax.scatter(Data['X'], Data['Y'], c='Blue', s=30.0, alpha=0.8, zorder=10)

	# Linear predict
	if ('LinRegLine' in Plot_Config):
		
		X_array = Data['X']
		Y_array, LinReg_Params, LinReg_pvalues, LinReg_Rsquareds = Calc_LinearPredict(Data['Y'], pd.DataFrame(Data['X']))
		ax.plot([np.min(X_array), np.max(X_array)], Y_array, color='Black', linewidth=2.8, alpha=0.6)

		Plot_Config['Text_1'] = 'Slope Estimator: {}\np-value={}, Rsq={}'.format(np.round(LinReg_Params.array[-1], 3), np.round(LinReg_pvalues.array[-1], 4), np.round(LinReg_Rsquareds, 3))

	return ax, New_PlotObject, Plot_Config

def Set_Configuration(ax, Data_Config, Plot_Config):

	XTicks = np.arange(Data_Config['Year_Start'], Data_Config['Year_End']+10, 10) \
									   if not ('XTicks' in Plot_Config) else Plot_Config['XTicks']
	XLabel = 'X'                       if not ('XLabel' in Plot_Config) else Plot_Config['XLabel']
	XLim   = [XTicks[0], XTicks[-1]]   if not ('XLim' in Plot_Config)   else Plot_Config['XLim']

	YTicks = np.arange(-50, 150, 50)   if not ('YTicks' in Plot_Config) else Plot_Config['YTicks']
	YLabel = 'Y'                       if not ('YLabel' in Plot_Config) else Plot_Config['YLabel']
	YLim   = [YTicks[0], YTicks[-1]]   if not ('YLim' in Plot_Config)   else Plot_Config['YLim']
	
	# X axis
	ax.set_xticks(XTicks)
	ax.tick_params(axis='x', labelsize=Plot_Config['Default_FontSize'])
	ax.set_xlabel(XLabel, fontsize=Plot_Config['Default_FontSize'])
	ax.set_xlim(XLim)

	# Y axis
	ax.set_yticks(YTicks)
	ax.tick_params(axis='y', labelsize=Plot_Config['Default_FontSize'])
	ax.set_ylabel(YLabel, fontsize=Plot_Config['Default_FontSize'])
	ax.set_ylim(YLim)

	ax.minorticks_on()
	ax.grid(which='major', color='grey', linewidth=0.6, alpha=0.6, zorder=15)
	ax.grid(which='minor', color='grey', linewidth=0.3, alpha=0.4, zorder=15)
	
	return ax

def Set_Text(ax, Plot_Config):

	LeftBottom_LineSpacing = 1.3
	LeftBottom_Position    = 0.02
	Text_FontSize          = 18

	if ('PlotTitle_Note' in Plot_Config):

		Figure_Title = 'VProfile, {VarLongName} {PlotTitle_Note}'.format(VarLongName=Plot_Config['VarLongName'], PlotTitle_Note=Plot_Config['PlotTitle_Note'])
		ax.set_title(Figure_Title, fontsize=10, pad=28)

	for i_Text in range(3):

		if ('Text_{}'.format(i_Text) in Plot_Config):
			
			ax.text(0.04, 1 - LeftBottom_Position * i_Text, r'{}'.format(Plot_Config['Text_' + str(i_Text)]), \
					fontsize=Text_FontSize, linespacing=LeftBottom_LineSpacing, horizontalalignment='left', verticalalignment='top', transform=ax.transAxes)
	
	return ax

def Save_Figure(ModelProperties, plt, Plot_Config, File_Config):

	if not ('FigurePath' in File_Config):
		
		FigurePath = 'Output_Figure/{ModelName}/ScatterPlot_{VarName}/'.format(ModelName=ModelProperties['ModelName'], VarName=Plot_Config['VarName'])
	
	else:

		FigurePath = File_Config['FigurePath']
	
	# Create directory if it's not existed
	if not os.path.exists(FigurePath):
		
		os.makedirs(FigurePath)
	
	FigurePath = FigurePath + 'ScatterPlot_{VarName}{FigureName_Note}.png'.format(VarName=Plot_Config['VarName'], FigureName_Note=Plot_Config['FigureName_Note'])
	
	# Save figure
	plt.tight_layout()
	plt.savefig(FigurePath, bbox_inches='tight')
	plt.close('all')

	return

def Plot(ModelProperties, Data, **kwargs):

	Data_Config           = kwargs.get('Data_Config', {})
	Plot_Config           = kwargs.get('Plot_Config', {})
	File_Config           = kwargs.get('File_Config', {})
	
	Plot_Config['Default_FontSize'] = 16
	
	# Create empty plot object dictionary
	PlotObject = {}

	# Plot
	fig, ax, Plot_Config = Figure_Initialization(Plot_Config)
	
	# Plot data
	ax, PlotObject['Plot_LinePlot'], Plot_Config = Draw_Plot(ax, Data, Data_Config, Plot_Config)

	# Plot config
	ax = Set_Configuration(ax, Data_Config, Plot_Config)

	# Plot text
	ax = Set_Text(ax, Plot_Config)

	# Save figure
	Save_Figure(ModelProperties, plt, Plot_Config, File_Config)

	return