# VProfile.py

# Import modules
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mplticker
import os

def Figure_Initialization(Plot_Config):

	# Create new figure
	FigureObject, AxesObject = plt.subplots(figsize=(8, 8), dpi=150)

	return FigureObject, AxesObject, Plot_Config

def Draw_Plot(ax, Data, Data_Config, Plot_Config, **kwargs):

	PlotMarker     = '.'             if not ('PlotMarker' in Plot_Config) else Plot_Config['PlotMarker']
	PlotMarkersize = 20              if not ('PlotMarkersize' in Plot_Config) else Plot_Config['PlotMarkersize']
	PlotLineWidth  = 6               if not ('PlotLineWidth' in Plot_Config) else Plot_Config['PlotLineWidth']
	
	PlotObject = []
	
	if isinstance(Data, dict):

		#for i_Line in range(len(Data['Data'])):
		for i_Line in [0, 2]:

			Legend_Labels = Plot_Config['Legend_Labels'][i_Line] if ('Legend_Labels' in Plot_Config) else None
			Data_Color    = Data['Data_Color'][i_Line]           if ('Data_Color' in Data) else None
			
			# Plot
			New_PlotObject = ax.plot(Data['Data'][i_Line][1:], Data_Config['Lev_Interp'][1:], \
									 color=Data_Color, label=Legend_Labels, alpha=0.7, \
									 marker=PlotMarker, markersize=PlotMarkersize, linewidth=PlotLineWidth, zorder=10)
			
			PlotObject.append(New_PlotObject)

			# Plot CI
			if ('Data_CI' in Data):

				ax.fill_betweenx(Data_Config['Lev_Interp'][1:], Data['Data_CI'][i_Line][0][1:], Data['Data_CI'][i_Line][-1][1:], \
								 fc=Data_Color, alpha=0.2, linewidth=PlotLineWidth, zorder=10)

	else:
		
		# Plot
		New_PlotObject = ax.plot(Data[1:], Data_Config['Lev_Interp'][1:], \
								 alpha=0.7, marker=PlotMarker, markersize=PlotMarkersize, linewidth=PlotLineWidth, zorder=10)
		
		PlotObject.append(New_PlotObject)
	
	if ('Data_Sig' in Data):

		for i_Sig in Data['Data_Sig']:

			ax.plot(Data['Data_Sig'][i_Sig][1:], Data_Config['Lev_Interp'][1:], \
					color='Grey', linestyle='None', alpha=0.9, \
					marker=i_Sig, markersize=10, zorder=15)
			PlotObject.append(New_PlotObject)

	if ('Legend_Labels' in Plot_Config):

		Legend = ax.legend(loc='upper right', prop={'size': 16})
		Legend.set_zorder(50)

	return ax, PlotObject

def Set_Configuration(ax, Plot_Config):

	XTicks = np.arange(-1000, 1500, 500)    if not ('XTicks' in Plot_Config) else Plot_Config['XTicks']
	XLabel = 'X axis'                       if not ('XLabel' in Plot_Config) else Plot_Config['XLabel']
	XLim   = \
	[np.min(XTicks), np.max(XTicks)]        if not ('XLim' in Plot_Config) else Plot_Config['XLim']
	XUnit  = ''                             if not ('XUnit' in Plot_Config) else r' ({})'.format(Plot_Config['XUnit'])

	YTicks = \
	np.concatenate((np.array([975]), np.arange(900, 0, -100))) if not ('YTicks' in Plot_Config) else Plot_Config['YTicks']
	YTicks_Minor = \
	np.arange(YTicks[0], YTicks[-1], 25)    if not ('YTicks_Minor' in Plot_Config) else Plot_Config['YTicks_Minor']
	YLabel = 'Pressure'                     if not ('YLabel' in Plot_Config) else Plot_Config['YLabel']
	YLim   = [975, 400]                     if not ('YLim' in Plot_Config) else Plot_Config['YLim']
	YUnit  = ' (hPa)'                       if not ('YUnit' in Plot_Config) else r' ({})'.format(Plot_Config['YUnit'])
	
	if ('HorizontalLines' in Plot_Config):

		for i_Line in range(len(Plot_Config['HorizontalLines'])):
			
			i_HorizontalLine = Plot_Config['HorizontalLines'][i_Line]
			ax.hlines(i_HorizontalLine['y'], *i_HorizontalLine['x'], linestyle='dashed')

	# X axis
	ax.set_xticks(XTicks)
	ax.tick_params(axis='x', labelsize=Plot_Config['Default_FontSize'])
	ax.set_xlabel(XLabel + XUnit, fontsize=Plot_Config['Default_FontSize'])
	ax.set_xlim(XLim)

	# Y axis
	ax.set_yscale('log')
	ax.set_yticks(YTicks)
	ax.set_yticks(YTicks_Minor, minor=True)
	ax.tick_params(axis='y', labelsize=Plot_Config['Default_FontSize'])
	ax.set_ylabel(YLabel + YUnit, fontsize=Plot_Config['Default_FontSize'])
	ax.set_ylim(YLim)

	ax.yaxis.set_major_formatter(mplticker.ScalarFormatter())
	ax.yaxis.get_major_formatter().set_scientific(False)
	ax.yaxis.get_major_formatter().set_useOffset(False)
	ax.yaxis.set_minor_formatter(mplticker.NullFormatter())
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
		
		FigurePath = 'Output_Figure/{ModelName}/VProfile_{VarName}/'.format(ModelName=ModelProperties['ModelName'], VarName=Plot_Config['VarName'])
	
	else:

		FigurePath = File_Config['FigurePath']
	
	# Create directory if it's not existed
	if not os.path.exists(FigurePath):
		
		os.makedirs(FigurePath)
	
	FigurePath = FigurePath + 'VProfile_{VarName}{FigureName_Note}.png'.format(VarName=Plot_Config['VarName'], FigureName_Note=Plot_Config['FigureName_Note'])
	
	# Save figure
	plt.tight_layout()
	plt.savefig(FigurePath, bbox_inches='tight')
	plt.close('all')

	return

def Plot(ModelProperties, Data, **kwargs):

	Data_Config           = kwargs.get('Data_Config', {})
	Plot_Config           = kwargs.get('Plot_Config', {})
	File_Config           = kwargs.get('File_Config', {})
	
	Plot_Config['Default_FontSize'] = 20
	
	# Create empty plot object dictionary
	PlotObject = {}

	# Plot
	fig, ax, Plot_Config = Figure_Initialization(Plot_Config)
	ax.vlines(0, 0, 1500, colors='Grey', linewidth=1.4, zorder=5)
	
	# Plot data
	ax, PlotObject['Plot_VProfile'] = Draw_Plot(ax, Data, Data_Config, Plot_Config)

	# Plot config
	ax = Set_Configuration(ax, Plot_Config)

	# Plot text
	ax = Set_Text(ax, Plot_Config)

	# Save figure
	Save_Figure(ModelProperties, plt, Plot_Config, File_Config)

	return