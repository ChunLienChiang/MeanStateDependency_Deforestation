# DataProcessing.py

import numpy as np
import pandas as pd
import netCDF4 as nc
import os

import scipy
import scipy.stats

g0 = 9.806

# Model properties

def Get_AbsolutePath(Path):

	if (Path[0]=='/'):

		AbsolutePath = Path

	else:

		AbsolutePath = os.path.abspath('../' + Path)

	return AbsolutePath

def Get_LandFraction(ModelProperties, **kwargs):

	Range        = kwargs.get('Range', 'All')

	Rawdata      = nc.Dataset(Get_AbsolutePath(ModelProperties['ReferenceFile_Path']['CTR']))

	if ('LANDFRAC' in Rawdata.variables):

		LandFraction = Rawdata.variables['LANDFRAC'][0, :, :]

	else:

		Rawdata      = nc.Dataset(Get_AbsolutePath('{ModelName}/Data/LANDFRAC/{ModelName}_0_LANDFRAC_CTR.nc'.format(ModelName=ModelProperties['ModelName'])))
		LandFraction = Rawdata.variables['LANDFRAC'][0, :, :]

	if (Range!='All'):
		
		data_lon_min, data_lon_max, data_lat_min, data_lat_max = Get_Range(ModelProperties, Range=Range, Type='Grid')
		LandFraction = LandFraction[data_lat_min:data_lat_max+1, data_lon_min:data_lon_max+1]

	return LandFraction

# ENSO

def Get_ONI_Data(ModelProperties, **kwargs):
	
	Data_TimeRange = kwargs.get('Data_TimeRange', ModelProperties['Properties']['Data_TimeRange'])
	Data_Frequency = kwargs.get('Data_Frequency', ModelProperties['Properties']['Data_Frequency'])
	
	# Obtain ONI data from file
	Data_ONI = pd.read_csv('ENSO_Data/ONI.csv', skiprows=0, header=0)
	Data_ONI = Data_ONI[(Data_ONI['Year'].between(Data_TimeRange[0], Data_TimeRange[1]+1))].to_numpy()[:, 1:].flatten()[6:-6]

	if (Data_Frequency=='Annually'):

		# Calculate annual mean
		Data_ONI = np.nanmean(Data_ONI.reshape(-1, 12), axis=1)
	
	return Data_ONI

def SplitData_ENSO(ModelProperties, Data, **kwargs):

	Data_TimeRange = kwargs.get('Data_TimeRange', ModelProperties['Properties']['Data_TimeRange'])
	ENSO_Intensity = kwargs.get('ENSO_Intensity', 'All')
	RM_Abnormal    = kwargs.get('RM_Abnormal', True)
	
	Data_SplitENSO = {}
	
	for i_ENSO in ['ElNino', 'LaNina', 'Neutral', 'All']:

		# Create tuple to store dimensions of new data
		NewAxis   = list(range(len(Data.shape)))
		NewAxis.pop(1)

		ENSO_Mask = Get_ENSO_Mask(Data_TimeRange, i_ENSO, data_f='Monthly', ENSO_Intensity=ENSO_Intensity, rm_atypical=RM_Abnormal)
		ENSO_mask = np.expand_dims(ENSO_Mask, axis=tuple(NewAxis))
		
		Data_SplitENSO[i_ENSO] = np.where(ENSO_mask, Data, np.nan)
	
	return Data_SplitENSO

def Get_ENSO_Mask(ModelProperties, ENSO_Phase, **kwargs):
	
	Data_TimeRange = kwargs.get('Data_TimeRange', ModelProperties['Properties']['Data_TimeRange'])
	Data_Frequency = kwargs.get('Data_Frequency', ModelProperties['Properties']['Data_Frequency']) # 'Monthly' or 'Annually'
	RM_Abnormal    = kwargs.get('RM_Abnormal', True)
	ENSO_Intensity = kwargs.get('ENSO_Intensity', 'All')
	ENSO_YearShift = kwargs.get('ENSO_YearShift', 0)
	ReturnType     = kwargs.get('ReturnType', None) # Set to 'Year' to return year number of ENSO events
	
	# Shift time range
	Data_TimeRange = list(np.array(Data_TimeRange)+ENSO_YearShift)
	
	# Return ONI Raw Data
	Data_ONI    = pd.read_csv('Data/ENSO/ONI.csv', skiprows=0, header=0)
	Data_ONI    = Data_ONI[(Data_ONI['Year'].between(Data_TimeRange[0], Data_TimeRange[1]))]
	Data_Events = Data_ONI.to_numpy()[:, 0]
	Data_Year   = Data_ONI.to_numpy()[:, 1]
	Data_ONI    = Data_ONI.to_numpy()[:, 2:].flatten()
	n_year      = Data_Events.shape[0]
	
	ENSO_mask = np.full((n_year), False)

	if (ENSO_Phase=='All'):

		ENSO_mask = np.full((n_year), True)

	elif ('Neutral' in ENSO_Phase):

		ENSO_mask[np.array([Events[-1] for Events in list(Data_Events)])=='N'] = True

		if (ENSO_Phase=='Neutral_1'):

			ENSO_mask_ONI = np.full((n_year), False)
			ENSO_mask_ONI[np.nanmean(Data_ONI.reshape((-1, 12)), axis=1)>0] = True
			ENSO_mask = ENSO_mask & ENSO_mask_ONI

		elif (ENSO_Phase=='Neutral_2'):

			ENSO_mask_ONI = np.full((n_year), False)
			ENSO_mask_ONI[np.nanmean(Data_ONI.reshape((-1, 12)), axis=1)<0] = True
			ENSO_mask = ENSO_mask & ENSO_mask_ONI

	else:

		if (ENSO_Intensity=='All'):

			ENSO_mask[np.array([Events[-1] for Events in list(Data_Events)])==ENSO_Phase[0]] = True

		else:

			for i_ENSO_Intensity in ENSO_Intensity:

				ENSO_mask[Data_Events=='{}{}'.format(i_ENSO_Intensity, ENSO_Phase[0])] = True
			
			if (RM_Abnormal==True) and (ENSO_Phase[0]=='E'):

				ENSO_mask[Data_Year==1987] = False
				ENSO_mask[Data_Year==1991] = False

			if (RM_Abnormal==True) and (ENSO_Phase[0]=='L'):

				ENSO_mask[Data_Year==1973] = False
				ENSO_mask[Data_Year==1975] = False
	
	if (ReturnType=='Year'):

		ENSO_mask = np.arange(Data_TimeRange[0], Data_TimeRange[1]+1)[ENSO_mask]
	
	if (Data_Frequency=='Monthly'):

		ENSO_mask = np.repeat(ENSO_mask, 12)
	
	return ENSO_mask

def Get_ENSO_MaskSet(ModelProperties, **kwargs):

	Data_TimeRange = kwargs.get('Data_TimeRange', ModelProperties['Properties']['Data_TimeRange'])
	Data_Frequency = kwargs.get('Data_Frequency', ModelProperties['Properties']['Data_Frequency'])
	RM_Abnormal    = kwargs.get('RM_Abnormal', True)

	ENSO_MaskSet = {}

	for i_ENSO in ['All', 'ElNino', 'LaNina', 'Neutral']:

		ENSO_MaskSet[i_ENSO] = list(Get_ENSO_Mask(ModelProperties, i_ENSO, Data_TimeRange=Data_TimeRange, Data_Frequency=Data_Frequency, \
												  RM_Abnormal=RM_Abnormal, ENSO_Intensity=['S', 'VS']))
	
	return ENSO_MaskSet

# Data mask

def Get_Land_Mask(ModelProperties, **kwargs):

	Range = kwargs.get('Range', 'MC_data_2')
	Type  = kwargs.get('Type')

	LandFraction = Get_LandFraction(ModelProperties, Range=Range)

	if (Type=='Ocean'):

		Land_Mask = (LandFraction==0)

	else:

		Land_Mask = (LandFraction>0.5)
	
	return Land_Mask

def Get_Offshore_Mask(ModelProperties, **kwargs):

	Range = kwargs.get('Range', 'MC_data_2')

	data_lon_min, data_lon_max, data_lat_min, data_lat_max = Get_Range(ModelProperties, Range=Range, Type='Grid')
	
	Offshore_Mask = np.loadtxt('../{}/Mask/Offshore.txt'.format(ModelProperties['ModelName']))
	Offshore_Mask = np.where(Offshore_Mask==1, True, False)
	
	return Offshore_Mask[data_lat_min:data_lat_max+1, data_lon_min:data_lon_max+1]

def Get_Range_Mask(ModelProperties, **kwargs):

	Range    = kwargs.get('Range', 'MC_data')
	Lon      = kwargs.get('Lon', Get_LonLat(ModelProperties, Range=Range)[0])
	Lat      = kwargs.get('Lat', Get_LonLat(ModelProperties, Range=Range)[1])

	data_lon_min, data_lon_max, data_lat_min, data_lat_max = Get_Range(ModelProperties, Range=Range, Type='Grid')

	Range_Mask = np.full((Lat.shape[0], Lon.shape[0]), False)
	Range_Mask[data_lat_min:data_lat_max+1, data_lon_min:data_lon_max+1] = True

	return Range_Mask

def Get_Year_Mask(ModelProperties, Mask_Year, **kwargs):

	Data_TimeRange = kwargs.get('Data_TimeRange', ModelProperties['Properties']['Data_TimeRange'])

	n_time = Data_TimeRange[1] - Data_TimeRange[0] + 1
	
	Year_array = np.full((n_time), False)
	Year_array[list(np.array(Mask_Year)-Data_TimeRange[0])] = True
	Year_array = np.repeat(Year_array, 12)

	return Year_array

def Get_LonLat(ModelProperties, **kwargs):

	ReferenceFile_Path = kwargs.get('ReferenceFile_Path', ModelProperties['ReferenceFile_Path']['CTR'])
	Range              = kwargs.get('Range', 'MC_data_2')

	Rawdata = nc.Dataset(Get_AbsolutePath(ReferenceFile_Path))
	
	Lon = Rawdata.variables['lon'][:]
	Lat = Rawdata.variables['lat'][:]
	
	if (Range!='All'):

		data_lon_min, data_lon_max, data_lat_min, data_lat_max = Get_Range(ModelProperties, Range=Range, Type='Grid')
		Lon = Lon[data_lon_min:data_lon_max+1]
		Lat = Lat[data_lat_min:data_lat_max+1]
	
	return Lon, Lat

def Get_Lev(ModelProperties):

	Rawdata      = nc.Dataset(Get_AbsolutePath(ModelProperties['ReferenceFile_Path']['CTR']))
	Lev          = Rawdata.variables['lev'][:]

	return Lev

def Get_Lev_Interp():

	return [1000, 975, 950, 925, 900, 875, 850, 825, 800, 775, 750, 725, 700, 650, 600, 550, 500, 400, 300, 200, 100]

def Get_Month_Mask():

	Month_Mask_Dict = {}
	Month_Name_List = ['J', 'A', 'S', 'O', 'N', 'D', 'J', 'F', 'M', 'A', 'M', 'J']
	
	for i_Month_length in range(2, 13):

		for i_Mask_index in range(13-i_Month_length):

			Month_Mask = [False] * 12

			for i_month in range(12):

				if (i_month in range(i_Mask_index, i_Mask_index+i_Month_length)):

					Month_Mask[i_month] = True

			Month_Mask_Dict[''.join([Month_Name_List[i] for i in range(i_Mask_index, i_Mask_index+i_Month_length)])] = Month_Mask
	
	Month_Mask_Dict['All']       = [True] * 12
	Month_Mask_Dict['July']      = [True, False, False, False, False, False, False, False, False, False, False, False]
	Month_Mask_Dict['August']    = [False, True, False, False, False, False, False, False, False, False, False, False]
	Month_Mask_Dict['September'] = [False, False, True, False, False, False, False, False, False, False, False, False]
	Month_Mask_Dict['October']   = [False, False, False, True, False, False, False, False, False, False, False, False]
	Month_Mask_Dict['November']  = [False, False, False, False, True, False, False, False, False, False, False, False]
	Month_Mask_Dict['December']  = [False, False, False, False, False, True, False, False, False, False, False, False]
	Month_Mask_Dict['January']   = [False, False, False, False, False, False, True, False, False, False, False, False]
	Month_Mask_Dict['February']  = [False, False, False, False, False, False, False, True, False, False, False, False]
	Month_Mask_Dict['March']     = [False, False, False, False, False, False, False, False, True, False, False, False]
	Month_Mask_Dict['April']     = [False, False, False, False, False, False, False, False, False, True, False, False]
	Month_Mask_Dict['May']       = [False, False, False, False, False, False, False, False, False, False, True, False]
	Month_Mask_Dict['June']      = [False, False, False, False, False, False, False, False, False, False, False, True]
	
	return Month_Mask_Dict

def Get_Lat_Weighting(ModelProperties, Range):
	
	Lon, Lat      = Get_LonLat(ModelProperties, Range=Range)
	LonLat_Weight = np.tile(np.cos(np.deg2rad(Lat))[:, None], Lon.shape[0])

	return LonLat_Weight

# ======================= CHOOSE PLOT RANGE =========================== #
def Get_Range(ModelProperties, **kwargs):

	Type = kwargs.get('Type')

	if (kwargs.get('Range')=='Nino3.4'):

		data_lon_min, data_lon_max, data_lat_min, data_lat_max = 190, 240, -5, 5

	elif (kwargs.get('Range')=='MC_plot'):

		data_lon_min, data_lon_max, data_lat_min, data_lat_max = 70, 180, -20, 20

	elif (kwargs.get('Range')=='MC_data'):

		data_lon_min, data_lon_max, data_lat_min, data_lat_max = 90, 150, -10, 10

	elif (kwargs.get('Range')=='MC_data_2'):
		
		# Default Maritime Continent domain (without east side of New Guinea)
		data_lon_min, data_lon_max, data_lat_min, data_lat_max = 90, 140, -10, 10

	elif (kwargs.get('Range')=='Borneo'):

		data_lon_min, data_lon_max, data_lat_min, data_lat_max = 108, 120, -5, 8

	elif (kwargs.get('Range')=='MC_TP1'):

		data_lon_min, data_lon_max, data_lat_min, data_lat_max = 70, 180, -5, 5

	elif (kwargs.get('Range')=='MC_TP2'):

		data_lon_min, data_lon_max, data_lat_min, data_lat_max = 70, 130, -5, 5

	elif (kwargs.get('Range')=='Tropical_Pacific_1'):

		data_lon_min, data_lon_max, data_lat_min, data_lat_max = 60, 300, -5, 5

	elif (kwargs.get('Range')=='Tropical_Pacific_2'):

		data_lon_min, data_lon_max, data_lat_min, data_lat_max = 60, 300, -40, 40

	elif (kwargs.get('Range')=='Tropical_Pacific_3'):

		# Default domain of variable preprocessing
		data_lon_min, data_lon_max, data_lat_min, data_lat_max = 50, 300, -40, 40

	elif (kwargs.get('Range')=='Test_1'):

		data_lon_min, data_lon_max, data_lat_min, data_lat_max = 100, 103, -6, -4

	elif (kwargs.get('Range')=='All'):

		data_lon_min, data_lon_max, data_lat_min, data_lat_max = 0, 360, -90, 90
	
	else:
		
		print('Main_Data_Range error')
		quit()
	
	if (Type=='Grid'):

		data_lon, data_lat = Get_LonLat(ModelProperties, Range='All')
		data_lon_min, data_lon_max, data_lat_min, data_lat_max = \
		np.argmin(abs(data_lon-data_lon_min)), np.argmin(abs(data_lon-data_lon_max)), np.argmin(abs(data_lat-data_lat_min)), np.argmin(abs(data_lat-data_lat_max))
	
	return data_lon_min, data_lon_max, data_lat_min, data_lat_max

# ========================== Calculation ============================== #

def Pressure_VI(scalar, plev):
	# plev must be arranged from small to large
	# Integrate along the first axis of scalar array
	return scipy.integrate.simps(scalar, plev, axis=0) / g0

# ======================== Data Processing ============================ #

def Calc_LonLat_Average(ModelProperties, Data, **kwargs):
	
	Range        = kwargs.get('Range', 'MC_data_2')
	Area         = kwargs.get('Area', 'Land')

	# Connect to lonlat_weight file
	LonLat_Weight = Get_Lat_Weighting(ModelProperties, Range)
	LandFraction  = Get_LandFraction(ModelProperties, Range=Range)
	
	if (Area=='Land'):

		Land_Mask = Get_Land_Mask(ModelProperties, Range=Range)

	elif (Area=='Ocean'):

		Land_Mask = Get_Land_Mask(ModelProperties, Range=Range, Type='Ocean')
		LandFraction = 1 - LandFraction

	elif (Area=='Offshore'):

		Land_Mask = Get_Offshore_Mask(ModelProperties, Range=Range)
		LandFraction = 1 - LandFraction
	
	# Combine masks
	Land_Mask = np.where(Land_Mask, 1, 0)
	Weighting = LonLat_Weight * LandFraction * Land_Mask

	# Calculate area weighted average
	if (len(Data.shape)==3):

		Weighting = np.tile(Weighting, (Data.shape[0], 1, 1))

	elif (len(Data.shape)==4):

		Weighting = np.tile(Weighting, (Data.shape[0], Data.shape[1], 1, 1))

	elif (len(Data.shape)==5):

		Weighting = np.tile(Weighting, (Data.shape[0], Data.shape[1], Data.shape[2], 1, 1))
	
	ma = np.ma.MaskedArray(Data, mask=np.isnan(Data))
	Weighting_Average = np.ma.average(ma, weights=Weighting, axis=(-2, -1))
	Weighting_Average = np.where(np.any(~np.isnan(Data), axis=(-2, -1)), Weighting_Average, np.nan)

	return Weighting_Average

def Calc_Month_Mask(Data, Month_Range, t_axis, **kwargs):

	Data_Frequency = kwargs.get('Data_Frequency', 'Monthly')

	if (Data_Frequency=='Monthly'):
		
		n_year = Data.shape[t_axis] // 12

	Month_Mask             = np.array(Get_Month_Mask()[Month_Range] * n_year)
	Dim_Month_Mask         = [None] * len(Data.shape)
	Dim_Month_Mask[t_axis] = slice(None)

	Data_Month_Mask        = np.where(Month_Mask[tuple(Dim_Month_Mask)], Data, np.nan)

	return Data_Month_Mask

def Calc_AnnualMean(Data, t_Axis):

	import warnings
	import skimage.measure

	with warnings.catch_warnings():

		warnings.simplefilter('ignore', category=RuntimeWarning)

		# Create new dimension tuple
		NewDim         = [1] * len(Data.shape)
		NewDim[t_Axis] = 12

		Data_AnnualMean = skimage.measure.block_reduce(Data, block_size=tuple(NewDim), func=np.nanmean, cval=np.nan)

	return Data_AnnualMean

def Calc_SeasonalCycle(Data, t_Axis, **kwargs):

	Return_Variability = kwargs.get('Return_Variability')
	Variability_Axis   = kwargs.get('Variability_Axis', (0, 1))

	# Calculate periodic cycle and tile along time dimension
	# Calculate dimension of climatology data
	Dim_Data         = list(Data.shape)
	Dim_Data         = Dim_Data[:t_Axis] + [Dim_Data[t_Axis]//12] + [12] + Dim_Data[t_Axis+1:]
	
	# Reshape the original data
	Data_Seasonal            = Data.reshape(tuple(Dim_Data))
	Data_Seasonal_Mean       = np.nanmean(Data_Seasonal, axis=t_Axis)
	_, Data_Seasonal_Mean_CI = Calc_CI(Data_Seasonal, axis=Variability_Axis)
	
	if (Return_Variability=='Std'):

		Data_Seasonal_Std = np.nanstd(Data_Seasonal, axis=Variability_Axis)
		
		return Data_Seasonal_Mean, Data_Seasonal_Mean_CI, Data_Seasonal_Std

	elif (Return_Variability=='Percentile'):

		Percentile_Pair = kwargs.get('Percentile_Pair')

		Data_Seasonal_Percentile = np.nanpercentile(Data_Seasonal, Percentile_Pair[0], axis=Variability_Axis) - \
								   np.nanpercentile(Data_Seasonal, Percentile_Pair[-1], axis=Variability_Axis)
		
		return Data_Seasonal_Mean, Data_Seasonal_Mean_CI, Data_Seasonal_Percentile
		
	else:

		return Data_Seasonal_Mean, Data_Seasonal_Mean_CI

def Calc_ClimAnomaly(Data, t_Axis, **kwargs):
	
	Calc_Cycle  = kwargs.get('Calc_Cycle', True)
	CyclePeriod = kwargs.get('CyclePeriod', 12)

	if (Calc_Cycle):
		
		# Calculate periodic cycle and tile along time dimension
		# Calculate dimension of climatology data
		Dim_Data_Clim    = list(Data.shape)
		Dim_Data_Clim    = Dim_Data_Clim[:t_Axis] + [Dim_Data_Clim[t_Axis]//CyclePeriod] + [CyclePeriod] + Dim_Data_Clim[t_Axis+1:]
		
		# Reshape the original data
		Dim_Tile         = [1] * len(Data.shape)
		Dim_Tile[t_Axis] = Data.shape[t_Axis] // CyclePeriod
		
		Data_Clim        = Data.reshape(tuple(Dim_Data_Clim))
		Data_Clim        = np.tile(np.nanmean(Data_Clim, axis=t_Axis).squeeze(), tuple(Dim_Tile))

	else:

		Data_Clim        = np.nanmean(Data, axis=t_Axis, keepdims=True)

	Data_ClimAnomaly = Data - Data_Clim

	return Data_ClimAnomaly

def Calc_Trend(Data, Time_Array, Time_axis):

	import statsmodels.api as sm

	Dim             = list(Data.shape)
	Dim[Time_axis]  = 1
	Dim             = tuple(Dim)

	Data_Trend = np.full(Dim, np.nan)
	Data_CI    = np.full(Dim, np.nan, dtype=object)

	for ind in np.ndindex(Dim):
		
		ind_t            = list(ind)
		ind_t[Time_axis] = slice(None)
		ind_t            = tuple(ind_t)
		
		LinReg = sm.OLS(Data[ind_t], sm.add_constant(Time_Array)).fit()

		Data_Trend[ind] = LinReg.params[-1]
		Data_CI[ind]    = LinReg.conf_int(alpha=0.05, cols=None)
	
	Data_Trend = Data_Trend.squeeze()
	Data_CI    = Data_CI.squeeze()

	return Data_Trend, Data_CI

def Calc_Detrend(Data, T_axis):
	
	global linear_trend

	Data          = np.swapaxes(Data, T_axis, 0)
	n_time        = Data.shape[0]
	x             = np.arange(n_time)
	
	Linear_Trend  = np.full_like(Data, np.nan)
	
	for i_ind, val in np.ndenumerate(Data[0, ...]):

		y   = Data[(slice(None),) + i_ind]
		idx = np.isfinite(x) & np.isfinite(y)
		Linear_Trend[(slice(None),) + i_ind] = np.polyval(np.polyfit(x[idx], y[idx], 1), x)
	
	Data_Detrend = np.swapaxes(Data - Linear_Trend, T_axis, 0)
	return Data_Detrend

def Calc_VI(Data, Lev, Integrate_Axis, **kwargs):

	import scipy.integrate

	Data    = np.where(np.isnan(Data), 0, Data)
	Data_VI = scipy.integrate.simps(Data, Lev, axis=Integrate_Axis) / g0

	return Data_VI

def Calc_CI(Data, **kwargs):

	Calc_Type = kwargs.get('Calc_Type', 'General')
	Calc_Axis = kwargs.get('Calc_Axis', (0, 1))

	if ((Calc_Type=='Diff') or (isinstance(Data, dict))):
	
		Data_1, Data_2 = Data['Data_1'], Data['Data_2']

		n_Data_1, Std_Data_1 = np.sum(~np.isnan(Data_1), axis=Calc_Axis), np.nanstd(Data_1, axis=Calc_Axis)
		n_Data_2, Std_Data_2 = np.sum(~np.isnan(Data_2), axis=Calc_Axis), np.nanstd(Data_2, axis=Calc_Axis)
		
		# Calculate the difference of two dataset
		Data_Mean = np.nanmean(Data_1, axis=Calc_Axis) - np.nanmean(Data_2, axis=Calc_Axis)

		# Calculate confidence interval
		sp      = np.sqrt((Std_Data_1**2)/(n_Data_1)+(Std_Data_2**2)/(n_Data_2))
		dof     = ((Std_Data_1**2)/(n_Data_1)+(Std_Data_2**2)/(n_Data_2))**2 / \
				  (((Std_Data_1**2)/(n_Data_1))**2/(n_Data_1-1)+((Std_Data_2**2)/(n_Data_2))**2/(n_Data_2-1))
		Data_CI = [Data_Mean - scipy.stats.t.ppf(1-0.025, dof)*sp, \
				   Data_Mean + scipy.stats.t.ppf(1-0.025, dof)*sp]
	
	elif (Calc_Type=='General'):
		
		n_Data, Std_Data = np.sum(~np.isnan(Data), axis=Calc_Axis), np.nanstd(Data, axis=Calc_Axis)

		# Calculate the difference of two dataset
		Data_Mean = np.nanmean(Data, axis=Calc_Axis)
		
		# Calculate confidence interval
		Data_CI = [Data_Mean - 1.96 * Std_Data/np.sqrt(n_Data), \
				   Data_Mean + 1.96 * Std_Data/np.sqrt(n_Data)]
	
	else:
		
		print('Argument \'calc_type\' should be \'Diff\' or \'General\'')
		return

	return Data_Mean, Data_CI

def Calc_pvalue(Data, **kwargs):

	import scipy.stats
	
	Calc_Type = kwargs.get('Calc_Type', 'General')
	Calc_Axis = kwargs.get('Calc_Axis', (0, 1))

	if ((Calc_Type=='Diff') or (isinstance(Data, dict))):
	
		Data_1, Data_2 = Data['Data_1'], Data['Data_2']
		
		if (len(Calc_Axis)>1):

			Dim_Flatten = np.delete(np.array(Data_1.shape), list(Calc_Axis))
			Dim_Flatten = np.concatenate((np.array([-1]), Dim_Flatten))
			
			Data_1, Data_2 = Data_1.reshape(tuple(Dim_Flatten)), Data_2.reshape(tuple(Dim_Flatten))
			Calc_Axis      = (0)

		# Calculate the difference of two dataset
		Data_Mean   = np.nanmean(Data_1, axis=Calc_Axis) - np.nanmean(Data_2, axis=Calc_Axis)
		Data_pvalue = scipy.stats.ttest_ind(Data_1, Data_2, axis=Calc_Axis, equal_var=False, nan_policy='omit').pvalue

	elif (Calc_Type=='General'):

		if (len(Calc_Axis)>1):

			Dim_Flatten = np.delete(np.array(Data.shape), list(Calc_Axis))
			Dim_Flatten = np.concatenate((np.array([-1]), Dim_Flatten))
			
			Data        = Data.reshape(tuple(Dim_Flatten))
			Calc_Axis   = (0)

		# Calculate the difference of two dataset
		Data_Mean   = np.nanmean(Data, axis=Calc_Axis)
		Data_pvalue = scipy.stats.ttest_1samp(Data, 0, axis=Calc_Axis, nan_policy='omit').pvalue
	
	else:
		
		print('Argument \'calc_type\' should be \'Diff\' or \'General\'')
		return

	return Data_Mean, Data_pvalue