# Map.py

# Import modules
import numpy as np
import matplotlib.pyplot as plt
import cartopy
import os

def Get_Extend(PlotColorbar):
	
	if (PlotColorbar['Levels'][0]==0):
		
		return 'max'

	elif (PlotColorbar['Levels'][-1]==0):

		return 'min'
		
	else:

		return 'both'

def Figure_Initialization(Data_Config):

	# Create new figure
	FigureObject, AxesObject = plt.subplots(figsize=(8, 5), dpi=200)
	plt.clf()
	
	# Set map range
	AxesObject               = plt.axes(projection=cartopy.crs.PlateCarree(central_longitude=180.0))
	AxesObject.set_extent([Data_Config['Lon'][0], Data_Config['Lon'][-1], Data_Config['Lat'][0], Data_Config['Lat'][-1]], crs=cartopy.crs.PlateCarree())
	
	return FigureObject, AxesObject

def Draw_Pcolormesh(ax, Data, Data_Config, Plot_Config):

	def Interp_Colormap(ColorList_Origin, n_ColorList_New):

		import scipy.interpolate as SciInterp

		n_ColorList_Origin = len(ColorList_Origin)
		ColorIndex_Origin  = np.linspace(0, 1, n_ColorList_Origin)
		ColorIndex_New     = np.linspace(0, 1, n_ColorList_New)

		Color_Interp       = SciInterp.interp1d(ColorIndex_Origin, np.array(ColorList_Origin), kind='linear', axis=0)
		ColorList_New      = Color_Interp(ColorIndex_New)

		return ColorList_New

	Plot_Alpha = 0.6 if not ('Plot_Alpha' in Plot_Config) else Plot_Config['Plot_Alpha']

	import matplotlib.colors
	import sys
	sys.path.append('..')
	import ClimPlot.Plot.NCLcmaps as NCLcmaps

	ColorList     = NCLcmaps.ColorList(Plot_Config['cmap'])
	ColorList_New = Interp_Colormap(ColorList, 17)
	cmap = matplotlib.colors.LinearSegmentedColormap.from_list(Plot_Config['cmap'], ColorList_New/256, N=len(ColorList_New))
	norm = matplotlib.colors.Normalize(vmin=Plot_Config['Colorbar']['Levels'][0]+0.5, vmax=Plot_Config['Colorbar']['Levels'][-1]-0.5)

	'''
	cmap = NCLcmaps.cmap(Plot_Config['cmap'])
	norm = matplotlib.colors.Normalize(vmin=Plot_Config['Colorbar']['Levels'][0], vmax=Plot_Config['Colorbar']['Levels'][-1])
	'''
	
	New_PlotObject = ax.pcolormesh(Data_Config['Lon'], Data_Config['Lat'], Data, shading='nearest', cmap=cmap, norm=norm, \
							 	   alpha=Plot_Alpha, zorder=10, vmin=Plot_Config['Colorbar']['Levels'][0]-0.5, vmax=Plot_Config['Colorbar']['Levels'][-1]+0.5, \
								   transform=cartopy.crs.PlateCarree())
	
	# Plot for SurfacePFT
	'''
	plt.colorbar(New_PlotObject, ticks=[Plot_Config['Colorbar']['Levels'][0], 0, Plot_Config['Colorbar']['Levels'][-1]], \
				 ax=ax, shrink=0.6, extend=Plot_Config['Colorbar']['Extend'], orientation='horizontal', pad=0.08)
	'''
	cbar = plt.colorbar(New_PlotObject, ticks=Plot_Config['Colorbar']['Levels'], ax=ax, shrink=0.6, orientation='vertical', pad=0.15)
	cbar.ax.tick_params(labelsize=8) 
	cbar.ax.set_yticklabels([r'1  - Desert, ice and ocean', \
							 r'2  - Needleleaf evergreen temperate tree', \
							 r'3  - Needleleaf evergreen boreal tree', \
							 r'4  - Needleleaf deciduous temperate tree', \
							 r'5  - Broadleaf evergreen tropical tree', \
							 r'6  - Broadleaf evergreen temperate tree', \
							 r'7  - Broadleaf deciduous tropical tree', \
							 r'8  - Broadleaf deciduous temperate tree', \
							 r'9  - Broadleaf deciduous boreal tree', \
							 r'10 - Broadleaf evergreen shrub', \
							 r'11 - Broadleaf deciduous temperate shrub', \
							 r'12 - Broadleaf deciduous boreal shrub', \
							 r'13 - C3 arctic grass', \
							 r'14 - C3 non-arctic grass', \
							 r'15 - C4 grass', \
							 r'16 - Corn', \
							 r'17 - Wheat'])

	return ax, New_PlotObject

def Set_Configuration(ax, Data_Config):

	import matplotlib.ticker as mt
	
	ax.coastlines(resolution='10m')
	ax.add_feature(cartopy.feature.NaturalEarthFeature('physical', 'land', '10m', edgecolor='face', facecolor='grey'))
	grid                = ax.gridlines(draw_labels=True, crs=cartopy.crs.PlateCarree(), linewidth=0.2, zorder=15)
	grid.top_labels     = False
	grid.xlabel_style   = {'size': 10}
	grid.ylabel_style   = {'size': 10}
	grid.xlocator       = mt.FixedLocator(np.mod(np.arange(Data_Config['Lon'][0], Data_Config['Lon'][-1]+30, 30)+180, 360)-180)
	
	return ax

def Set_Text(ax, Plot_Config):

	LeftBottom_LineSpacing = 1.3
	LeftBottom_Position    = 0.025
	Text_FontSize          = 10

	if ('PlotTitle_Note' in Plot_Config):

		Figure_Title = 'Map, {VarLongName} {PlotTitle_Note}'.format(VarLongName=Plot_Config['VarLongName'], PlotTitle_Note=Plot_Config['PlotTitle_Note'])
		ax.set_title(Figure_Title, fontsize=10, pad=28)

	for i_Text in range(3):

		if ('Text_{}'.format(i_Text) in Plot_Config):
			
			ax.text(0, 1 + LeftBottom_Position * i_Text, r'{}'.format(Plot_Config['Text_' + str(i_Text)]), \
					fontsize=Text_FontSize, linespacing=LeftBottom_LineSpacing, horizontalalignment='left', verticalalignment='bottom', transform=ax.transAxes)
	
	return ax

def Save_Figure(ModelProperties, plt, Plot_Config, File_Config):

	if not ('FigurePath' in File_Config):
		
		FigurePath = 'Output_Figure/{ModelName}/Map_{VarName}/'.format(ModelName=ModelProperties['ModelName'], VarName=Plot_Config['VarName'])
	
	else:

		FigurePath = File_Config['FigurePath']
	
	# Create directory if it's not existed
	if not os.path.exists(FigurePath):
		
		os.makedirs(FigurePath)
	
	FigurePath = FigurePath + 'Map_{VarName}{FigureName_Note}.png'.format(VarName=Plot_Config['VarName'], FigureName_Note=Plot_Config['FigureName_Note'])
	
	# Save figure
	plt.tight_layout()
	plt.savefig(FigurePath, bbox_inches='tight')
	plt.close('all')

	return

def Plot(ModelProperties, Data, **kwargs):

	Data_Config           = kwargs.get('Data_Config', {})
	Plot_Config           = kwargs.get('Plot_Config', {})
	File_Config           = kwargs.get('File_Config', {})
	
	# Get colorbar extend
	Plot_Config['Colorbar']['Extend'] = Get_Extend(Plot_Config['Colorbar'])
	
	# Create empty plot object dictionary
	PlotObject = {}

	# Plot
	fig, ax = Figure_Initialization(Data_Config)

	# Plot data
	ax, PlotObject['Pcolormesh_SST'] = Draw_Pcolormesh(ax, Data, Data_Config, Plot_Config)

	# Plot config
	ax = Set_Configuration(ax, Data_Config)

	# Plot text
	ax = Set_Text(ax, Plot_Config)

	# Save figure
	Save_Figure(ModelProperties, plt, Plot_Config, File_Config)

	return