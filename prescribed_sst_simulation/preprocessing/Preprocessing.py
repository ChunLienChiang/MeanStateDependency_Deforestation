"""
tool_ClimData_Preprocessing.py
==========================
The tools used to process the climate dataset.
"""

import numpy as np
import xarray as xr
import scipy.integrate
import os
import json

# Declare constants
Cp = 1006
g  = 9.8
Lv = 2.25e6

# Get configuration
Config = json.load(open(os.path.join(os.path.dirname(__file__), '../config.json')))

def Get_Lev_Interp():

	"""
	Get the level array to be interpolated
	"""

	return np.array([1000, 975, 950, 925, 900, 875, 850, 825, 800, 775, 750, 725, 700, 650, 600, 550, 500, 400, 300, 200, 100, 90, 80, 70, 60, 50])

def Get_RefData(\
		Component='atm', \
		Var=None, \
	):
	
	"""
	Connect to the reference dataset. The path of required dataset is provided in config.json
	==========================
	Argument:

		Component (str): the component model to get reference file
		
		Var (str): the variable to be read

	Output:

		Data (numpy array): the data in the nc file
	==========================
	"""
	
	# Check whether argument Var is empty
	if (Var is None):

		raise ValueError('The argument Var must not be empty.')

	# Get reference file
	Ref_File = Config['Ref_File'][Component]

	# Connect to nc file
	Data = xr.open_dataset(Ref_File)[Var].to_numpy()

	return Data

def Get_RefInfo():
	
	"""
	Get the reference information from the reference dataset. The path of required reference dataset is provided in config.json
	==========================
	Output:

		RefInfo (dict): the reference information
	==========================
	"""

	RefInfo = {\
		'Lat': Get_RefData(Var='lat'), \
		'Lon': Get_RefData(Var='lon'), \
		'LandFraction': Get_RefData(Var='LANDFRAC'), \
	}
	
	return RefInfo

def Crop_Range(Data, Range, Range_Original='Global_Analysis'):

	"""
	Crop the data to fit the given range
	==========================
	Argument:

		Data (numpy array): the data to be cropped

		Range (str): the range to be fitted

		Range_Original (str): optional, the original range of data. Default is global

	Output:

		Data_Crop (tuple): cropped data
	==========================
	"""

	Lon = Get_RefData(Var='lon')
	Lat = Get_RefData(Var='lat')

	if (tuple(Data.shape[-2:]) != tuple([len(i) for i in Crop_Lat_Lon(Range_Original)])):

		raise ValueError('Error in Crop_Range: given array does not meet original range.')
	
	else:
		
		Lat_Min, Lat_Max, Lon_Min, Lon_Max = Get_RangeBoundary(Range)
		Arg_Lon_Min = np.argmin(np.abs(Lon-Lon_Min))
		Arg_Lon_Max = np.argmin(np.abs(Lon-Lon_Max)) + 1
		Arg_Lat_Min = np.argmin(np.abs(Lat-Lat_Min))
		Arg_Lat_Max = np.argmin(np.abs(Lat-Lat_Max)) + 1
	
		if (Range_Original != 'Global_Analysis'):

			Lon_Min_Original, Lon_Max_Original, Lat_Min_Original, Lat_Max_Original = Get_RangeBoundary(Range_Original)
			Arg_Lon_Min_Original = np.argmin(np.abs(Lon-Lon_Min_Original))
			Arg_Lon_Max_Original = np.argmin(np.abs(Lon-Lon_Max_Original)) + 1
			Arg_Lat_Min_Original = np.argmin(np.abs(Lat-Lat_Min_Original))
			Arg_Lat_Max_Original = np.argmin(np.abs(Lat-Lat_Max_Original)) + 1

			Arg_Lon_Min          = Arg_Lon_Min - Arg_Lon_Min_Original
			Arg_Lon_Max          = Arg_Lon_Max - Arg_Lon_Max_Original
			Arg_Lat_Min          = Arg_Lat_Min - Arg_Lat_Min_Original
			Arg_Lat_Max          = Arg_Lat_Max - Arg_Lat_Max_Original
		
		if ((Arg_Lat_Max in Lat[[0, -1]]) and (Arg_Lon_Max in Lon[[0, -1]])):

			Data_Crop = Data[..., Arg_Lat_Min:, Arg_Lon_Min:]
		
		elif (Arg_Lat_Max in Lat[[0, -1]]) and ~(Arg_Lon_Max in Lon[[0, -1]]):

			Data_Crop = Data[..., Arg_Lat_Min:, Arg_Lon_Min:Arg_Lon_Max]
		
		elif ~(Arg_Lat_Max in Lat[[0, -1]]) and (Arg_Lon_Max in Lon[[0, -1]]):

			Data_Crop = Data[..., Arg_Lat_Min:Arg_Lat_Max, Arg_Lon_Min:]
		
		else:

			Data_Crop = Data[..., Arg_Lat_Min:Arg_Lat_Max, Arg_Lon_Min:Arg_Lon_Max]
		
		return Data_Crop

def Crop_Lat_Lon(Range):

	"""
	Crop the latitude and longitude data to fit the given range
	==========================
	Argument:

		Range (str): the range to be fitted

	Output:

		Lat_Crop (tuple): cropped latitude

		Lon_Crop (tuple): cropped longitude
	==========================
	"""
	
	Lon = Get_RefData(Var='lon')
	Lat = Get_RefData(Var='lat')

	Lat_Min, Lat_Max, Lon_Min, Lon_Max = Get_RangeBoundary(Range)
	Arg_Lon_Min = np.argmin(np.abs(Lon-Lon_Min))
	Arg_Lon_Max = np.argmin(np.abs(Lon-Lon_Max))
	Arg_Lat_Min = np.argmin(np.abs(Lat-Lat_Min))
	Arg_Lat_Max = np.argmin(np.abs(Lat-Lat_Max))

	Lat_Crop = Lat[Arg_Lat_Min:Arg_Lat_Max+1]
	Lon_Crop = Lon[Arg_Lon_Min:Arg_Lon_Max+1]

	return Lat_Crop, Lon_Crop

def Get_RangeBoundary(Range):

	"""
	Get the range boundary
	==========================
	Argument:

		Range (str): the range to get boundary

	Output:

		RangeBoundary (tuple): (minimum latitude, maximum latitude, minimum longitude, maximum longitude), respectively
	==========================
	"""

	if (Range == 'MC_Analysis'):
		
		return -10, 10, 90, 140
	
	elif (Range == 'MC_Extended_Analysis'):

		return -15, 15, 85, 145

	elif (Range == 'Borneo_Analysis'):

		return -5, 8, 108, 120

	elif (Range == 'Global_Analysis'):

		return -90, 90, 0, 360

	elif (Range == 'NH_Analysis'):

		return 0, 90, 0, 360

	elif (Range == 'SH_Analysis'):

		return -90, 0, 0, 360

	else:
		
		raise ValueError('The argument Range is not available.')

def Get_RangeMask(Range, Lat=None, Lon=None):

	"""
	Get the range mask
	==========================
	Argument:

		Range (str): the range to generate range mask

		Lat (numpy array): latitude. If there's no value provided, the latitude in the reference dataset would be used

		Lon (numpy array): longitude. If there's no value provided, the longitude in the reference dataset would be used
	
	Output:

		RangeMask (numpy array): range mask
	==========================
	"""

	# Check whether the arguments Lat and Lon are given
	if ((Lat is None) or (Lon is None)):
		
		RefInfo = Get_RefInfo()
		Lat     = RefInfo['Lat']
		Lon     = RefInfo['Lon']
	
	# Get the range boundary
	Lat_Min, Lat_Max, Lon_Min, Lon_Max = Get_RangeBoundary(Range)

	# Create the range mask
	RangeMask = np.where((Lat[:, None]>=Lat_Min)&(Lat[:, None]<=Lat_Max)&(Lon[None, :]>=Lon_Min)&(Lon[None, :]<=Lon_Max), True, False)

	return RangeMask

def Get_LandMask(LandFraction=None, MaskType='Land'):
	
	"""
	Get land mask from land fraction data.
	==========================
	Argument:

		LandFraction (numpy array): the land fraction data. If there's no value provided, the land fraction in the reference dataset would be used

		MaskType (str): 'Land' for land region (default); 'Ocean' for ocean region

	Output:

		LandMask (numpy array): the mask of land (or ocean)
	==========================
	"""
	
	# Check if the argument LandFraction is empty
	if (LandFraction is None):

		LandFraction = Get_RefInfo()['LandFraction']

	# Determine what tpye of mask should be return
	if (MaskType == 'Land'):

		LandMask = np.where(LandFraction>=0.5, True, False)

	elif (MaskType == 'Ocean'):

		LandMask = np.where(LandFraction<0.5, True, False)
	
	else:

		raise ValueError('The argument MaskType must be Land or Ocean')

	return LandMask

def Get_SeasonMask(Season, Data_Shift=True):

	"""
	Get the monthly mask of the given season
	==========================
	Argument:

		Season (str): the season to get mask. 'DJF', 'JA' are available

		Data_Shift (bool): whether the mask should be shifted to match data (begin at July). Default is True

	Output:

		SeasonMask (numpy array): the mask of the given season
	==========================
	"""
	
	if (Season == 'DJF'):

		SeasonMask = np.array([True, True, False, False, False, False, False, False, False, False, False, True])

	elif (Season == 'JA'):

		SeasonMask = np.array([False, False, False, False, False, False, True, True, False, False, False, False])

	else:

		raise ValueError('The argument Season is not available.')

	if (Data_Shift):
		
		SeasonMask = np.roll(SeasonMask, 6)
	
	return SeasonMask

def Calc_SpatialAverage(\
		Data, \
		Spatial_Axis=(-2, -1), \
		Bool_LatWeighted=True, \
		LandMask='Land', \
		RangeMask=None, \
		Range_Original='Global_Analysis', \
	):

	"""
	Calculate the spatial average of the data
	==========================
	Argument:

		Data (numpy array)

		Spatial_Axis (tuple of int): optional. The axis numbers of spatial (latitude and longtitude) dimension. Default is (-2, -1) (the last two dimensions)

		Bool_LatWeighted (boolean): Default is True. If the argument is True, the spatial average would be calculated considering latitude weighted

		LandMask (str): Default is None. If the argument is given, the land mask would be used to filter data

		RangeMask (str): Default is None. If the argument is given, the range mask would be used to filter data

		Range_Original (str): optional, the original range of data. Default is global
	
	Output:

		SpatialAverage (numpy array)
	==========================
	"""

	# Check whether the dimensions are correct
	if not (isinstance(Spatial_Axis, tuple)):

		raise ValueError('The argument Spatial_Axis must be a tuple.')
	
	if (len(Spatial_Axis) != 2):

		raise ValueError('The length of the argument Spatial_Axis must be 2.')

	# Create latitude weights
	Weight_Lat = np.cos(np.deg2rad(Get_RefInfo()['Lat']))[:, None]
	
	# Create land mask weights
	if (LandMask is None):
		
		Weight_LandMask = np.full((Data.shape[Spatial_Axis[0]], Data.shape[Spatial_Axis[1]]), 1)

	else:

		Weight_LandMask = np.where(Get_RefInfo()['LandFraction']>=0.5, Get_RefInfo()['LandFraction'], 0)
	
	# Create range mask weights
	if (RangeMask is None):
		
		Weight_RangeMask = np.full((Data.shape[Spatial_Axis[0]], Data.shape[Spatial_Axis[1]]), 1)

	else:

		Weight_RangeMask = np.where(Get_RangeMask(RangeMask), 1, 0)
	
	# Calculate spatial average
	if (Bool_LatWeighted):

		Weights = Crop_Range(Weight_Lat*Weight_LandMask*Weight_RangeMask, Range_Original)
	
	else:

		Weights = Crop_Range(Weight_LandMask*Weight_RangeMask, Range_Original)
	
	if ((np.ndim(Weights) - np.ndim(Data) == 1)):

		Weights = Weights[None, ...]
	
	Weights = np.broadcast_to(Weights, Data.shape)
	
	SpatialAverage = np.ma.average(np.ma.MaskedArray(Data, mask=np.isnan(Data)), axis=Spatial_Axis, weights=Weights)
		
	return SpatialAverage

def Calc_SeasonalCycle(Data, Time_Axis=0, Data_Frequency='Monthly'):

	"""
	Calculate the seasonal cycle of the data
	==========================
	Argument:

		Data (numpy array)

		Time_Axis (int): optional. The axis number of time dimension. Default is 0 (the first dimension)

		Data_Frequency (str): optional. 'Monthly' or 'Daily'. The frequency of time dimension of the data
	
	Output:

		SeasonalCycle (numpy array)
	==========================
	"""

	# Check whether the dimensions are correct
	if (not isinstance(Time_Axis, int)):

		raise ValueError('The argument Time_Axis must not integer.')
	
	# Set the period of the seasonal cycle
	if (Data_Frequency == 'Monthly'):

		Time_Period = 12

	if (Data_Frequency == 'Daily'):

		Time_Period = 365
	
	# Check whether the time dimension is a multiple of 12 (or 365) if Data_Frequency is monthly (or daily)
	if (Data.shape[Time_Axis] % Time_Period != 0) & (Data_Frequency == 'Monthly'):

		raise ValueError(\
			'If the argument Data_Frequency is monthly, '\
			'the time dimension of Data must be a multiple of 12. '\
			'If the argument Data_Frequency is daily, '\
			'the time dimension of Data must be a multiple of 365.'\
		)

	# Calculate seasonal cycle
	SeasonalCycle = np.nanmean(Data.reshape(tuple([*Data.shape[:Time_Axis], Data.shape[Time_Axis]//Time_Period, Time_Period, *Data.shape[Time_Axis+1:]])), axis=Time_Axis)

	return SeasonalCycle

def Calc_AnnualMean(Data, Time_Axis=0, Data_Frequency='Monthly', Season=None):

	"""
	Calculate the annual mean of the data
	==========================
	Argument:

		Data (numpy array)

		Time_Axis (int): optional. The axis number of time dimension. Default is 0 (the first dimension)

		Data_Frequency (str): optional. 'Monthly' or 'Daily'. The frequency of time dimension of the data
	
		Season (str): optional. The season to mask. Default is None (no mask). 'DJF', 'JA' are available
	
	Output:

		AnnualMean (numpy array)
	==========================
	"""

	# Check whether the dimensions are correct
	if (not isinstance(Time_Axis, int)):

		raise ValueError('The argument Time_Axis must not integer.')
	
	# Set the period of the annual mean
	if (Data_Frequency == 'Monthly'):

		Time_Period = 12

	if (Data_Frequency == 'Daily'):

		Time_Period = 365
	
	# Check whether the time dimension is a multiple of 12 (or 365) if Data_Frequency is monthly (or daily)
	if (Data.shape[Time_Axis] % Time_Period != 0) & (Data_Frequency == 'Monthly'):

		raise ValueError(\
			'If the argument Data_Frequency is monthly, '\
			'the time dimension of Data must be a multiple of 12. '\
			'If the argument Data_Frequency is daily, '\
			'the time dimension of Data must be a multiple of 365.'\
		)
	
	# Mask the season
	if (Season is not None):

		# Get season mask. Tile the season mask to the shape of Data
		dim_tile               = [1] * np.ndim(Data)
		dim_tile[Time_Axis]    = Data.shape[Time_Axis]//Time_Period
		dim_reshape            = [1] * np.ndim(Data)
		dim_reshape[Time_Axis] = 12
		Data = np.where(np.tile(Get_SeasonMask(Season).reshape(tuple(dim_reshape)), dim_tile), Data, np.nan)
	
	# Calculate annual mean
	AnnualMean = np.nanmean(Data.reshape(tuple([*Data.shape[:Time_Axis], Data.shape[Time_Axis]//Time_Period, Time_Period, *Data.shape[Time_Axis+1:]])), axis=Time_Axis+1)

	return AnnualMean

def Calc_VI(Data, Lev, Integrate_Axis):

	"""
	Calculate the vertical integration
	==========================
	Argument:

		Data (numpy array)

		Lev (numpy array): vertical level

		Integrate_Axis (int): the axis number of vertical level dimension

	Output:

		VI (numpy array)
	==========================
	"""

	Data = np.where(np.isnan(Data), 0, Data)
	VI   = scipy.integrate.simps(Data, Lev, axis=Integrate_Axis) / g

	return VI

def Calc_CI(Data, **kwargs):

	"""
	Calculate the confidence interval
	==========================
	Argument:

		Data (numpy array)

		Data_Type (numpy array): (optional) the type of t-test. 'General' for difference from zero; 'Diff' for difference of two group means. Default is 'General'
			If the argument if 'Diff', the data must be a list that contain two numpy arrays

		Test_Axis (int): (optional) the axis number of vertical level dimension. Default is (0, 1)
		
		Bool_ReturnWidth (boolean): (optional) whether return the width of CI

	Output:

		VI (numpy array)
	==========================
	"""

	Data_Type        = kwargs.get('Data_Type', 'General')
	Test_Axis        = kwargs.get('Test_Axis', (0, 1))
	Bool_ReturnWidth = kwargs.get('Bool_ReturnWidth', False)

	if ((Data_Type=='Diff') or (isinstance(Data, dict))):
	
		Data_1, Data_2 = Data['Data_1'], Data['Data_2']

		n_Data_1, Std_Data_1 = np.sum(~np.isnan(Data_1), axis=Test_Axis), np.nanstd(Data_1, axis=Test_Axis)
		n_Data_2, Std_Data_2 = np.sum(~np.isnan(Data_2), axis=Test_Axis), np.nanstd(Data_2, axis=Test_Axis)
		
		# Calculate the difference of two dataset
		Data_Mean = np.nanmean(Data_1, axis=Test_Axis) - np.nanmean(Data_2, axis=Test_Axis)

		# Calculate confidence interval
		sp      = np.sqrt((Std_Data_1**2)/(n_Data_1)+(Std_Data_2**2)/(n_Data_2))
		dof     = ((Std_Data_1**2)/(n_Data_1)+(Std_Data_2**2)/(n_Data_2))**2 / \
			(((Std_Data_1**2)/(n_Data_1))**2/(n_Data_1-1)+((Std_Data_2**2)/(n_Data_2))**2/(n_Data_2-1))
		Data_CI = [\
			Data_Mean - scipy.stats.t.ppf(1-0.025, dof)*sp, \
			Data_Mean + scipy.stats.t.ppf(1-0.025, dof)*sp, \
		]
		Data_CIWidth = scipy.stats.t.ppf(1-0.025, dof)*sp
	
	elif (Data_Type=='General'):
		
		n_Data, Std_Data = np.sum(~np.isnan(Data), axis=Test_Axis), np.nanstd(Data, axis=Test_Axis)

		# Calculate the difference of two dataset
		Data_Mean = np.nanmean(Data, axis=Test_Axis)
		
		# Calculate confidence interval
		Data_CI = [\
			Data_Mean - 1.96 * Std_Data/np.sqrt(n_Data), \
			Data_Mean + 1.96 * Std_Data/np.sqrt(n_Data), \
		]
		Data_CIWidth = 1.96 * Std_Data/np.sqrt(n_Data)

	if (Bool_ReturnWidth):

		return Data_Mean, Data_CI, Data_CIWidth
	
	else:

		return Data_Mean, Data_CI