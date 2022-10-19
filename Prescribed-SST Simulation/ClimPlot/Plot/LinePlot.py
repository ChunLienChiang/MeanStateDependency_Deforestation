# LinePlot.py

# Import modules
import numpy as np
import matplotlib.pyplot as plt
import os


def Figure_Initialization(Plot_Config):

	# Create new figure
	FigureObject, AxesObject = plt.subplots(figsize=(8, 8), dpi=150)

	return FigureObject, AxesObject, Plot_Config

def Draw_Plot(ax, Data, Data_Config, Plot_Config, **kwargs):

	Plot_XArray    = np.arange(Data_Config['Year_Start'], Data_Config['Year_End']+1) \
									 if not ('XArray' in Plot_Config) else Plot_Config['XArray']
	PlotMarker     = '.'             if not ('PlotMarker' in Plot_Config) else Plot_Config['PlotMarker']
	PlotMarkersize = 20              if not ('PlotMarkersize' in Plot_Config) else Plot_Config['PlotMarkersize']
	PlotLineWidth  = 6               if not ('PlotLineWidth' in Plot_Config) else Plot_Config['PlotLineWidth']
	
	PlotColor      = kwargs.get('PlotColor', None)

	# Plot
	if not isinstance(Data, dict):

		print('Data argument must be in type of dictionary')
		
		return

	else:

		if ('LinePlot' in Data):

			for i_Line in range(len(Data['LinePlot'])):

				if (PlotColor is None):

					New_PlotObject = ax.plot(Plot_XArray, Data['LinePlot'][i_Line], alpha=0.7, \
											marker=PlotMarker, markersize=PlotMarkersize, linewidth=PlotLineWidth, zorder=10)
				
				else:

					New_PlotObject = ax.plot(Plot_XArray, Data['LinePlot'][i_Line], alpha=0.7, color=PlotColor[i_Line], \
											marker=PlotMarker, markersize=PlotMarkersize, linewidth=PlotLineWidth, zorder=10)

	return ax, New_PlotObject

def Set_Configuration(ax, Data_Config, Plot_Config):

	XTicks          = np.arange(Data_Config['Year_Start'], Data_Config['Year_End']+10, 10) \
			                					if not ('XTicks' in Plot_Config) else Plot_Config['XTicks']
	XLabel          = 'X axis'                  if not ('XLabel' in Plot_Config) else Plot_Config['XLabel']
	XLim            = [Data_Config['Year_Start'], Data_Config['Year_End']] \
						        				if not ('XLim' in Plot_Config) else Plot_Config['XLim']

	YTicks          = np.arange(-2.5, 3.0, 0.5) if not ('YTicks' in Plot_Config) else Plot_Config['YTicks']
	YLabel          = 'Y axis'                  if not ('YLabel' in Plot_Config) else Plot_Config['YLabel']
	YLim            = [-2.5, 2.5]               if not ('YTicks' in Plot_Config) else Plot_Config['YLim']
	
	HorizontalLines = []                        if not ('HorizontalLines' in Plot_Config) else Plot_Config['HorizontalLines']
	
	# Zero lines
	ax.hlines(0, *XLim, color='Grey')
	ax.vlines(0, *YLim, color='Grey')

	for i_Line in HorizontalLines:

		ax.hlines(i_Line, *XLim, color='Grey', linestyle='dashed')

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
	Text_FontSize          = 10

	if ('PlotTitle_Note' in Plot_Config):

		Figure_Title = 'VProfile, {VarLongName} {PlotTitle_Note}'.format(VarLongName=Plot_Config['VarLongName'], PlotTitle_Note=Plot_Config['PlotTitle_Note'])
		ax.set_title(Figure_Title, fontsize=10, pad=28)

	for i_Text in range(3):

		if ('Text_{}'.format(i_Text) in Plot_Config):
			
			ax.text(0, 1 + LeftBottom_Position * i_Text, r'{}'.format(Plot_Config['Text_' + str(i_Text)]), \
					fontsize=Text_FontSize, linespacing=LeftBottom_LineSpacing, horizontalalignment='left', verticalalignment='bottom', transform=ax.transAxes)
	
	return ax

def Save_Figure(ModelProperties, plt, Plot_Config, File_Config):

	if not ('FigurePath' in File_Config):
		
		FigurePath = 'Output_Figure/{ModelName}/LinePlot_{VarName}/'.format(ModelName=ModelProperties['ModelName'], VarName=Plot_Config['VarName'])
	
	else:

		FigurePath = File_Config['FigurePath']
	
	# Create directory if it's not existed
	if not os.path.exists(FigurePath):
		
		os.makedirs(FigurePath)
	
	FigurePath = FigurePath + 'LinePlot_{VarName}{FigureName_Note}.png'.format(VarName=Plot_Config['VarName'], FigureName_Note=Plot_Config['FigureName_Note'])
	
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
	ax.vlines(0, 0, 1500, colors='Grey', linewidth=1.4, zorder=5)
	
	# Plot data
	ax, PlotObject['Plot_LinePlot'] = Draw_Plot(ax, Data, Data_Config, Plot_Config)
	
	# Plot config
	ax = Set_Configuration(ax, Data_Config, Plot_Config)

	# Plot text
	ax = Set_Text(ax, Plot_Config)

	# Save figure
	Save_Figure(ModelProperties, plt, Plot_Config, File_Config)

	return