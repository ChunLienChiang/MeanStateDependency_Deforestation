# DataProcessing.py

def Get_Model_Properties():

	import numpy as np
	import netCDF4 as nc

	Ref_Data = nc.Dataset('/work/data/CMIP6_B_piControl/TS/b.e21.B1850.f09_g17.CMIP6-piControl.001.cam.h0.TS.000101-009912.nc')

	Model_Properties = {}

	Model_Properties['Lat'] = Ref_Data.variables['lat'][:]
	Model_Properties['Lon'] = Ref_Data.variables['lon'][:]
	Model_Properties['Lev'] = Ref_Data.variables['lev'][:]
	Model_Properties['Lev_Interp']   = np.array([1000, 975, 950, 925, 900, 875, 850, 825, 800, 775, 750, 725, 700, 650, 600, 550, 500, 400, 300, 200, 100])

	Model_Properties['LandFraction'] = np.loadtxt('Data/LANDFRAC/LANDFRAC.txt')

	return Model_Properties

def Get_Range(Range):

	if (Range == 'Global_Analysis'):
		
		return 0, 360, -90, 90

	elif (Range == 'MC_Analysis'):
		
		return 90, 140, -10, 10

	elif (Range == 'MC_Analysis_Extended'):
		
		return 85, 145, -15, 15
	
	elif (Range == 'Borneo_Analysis'):
		
		return 108, 120, -5, 7.5

	elif (Range == 'Nino34_Analysis'):
		
		return 190, 240, -5, 5

	elif (Range == 'PDO_Analysis'):
		
		return 120, 250, 20, 90

	elif (Range == 'PDO_NorthPacific_Analysis'):
		
		return 140, 200, 30, 50

	elif (Range == 'AMO_Analysis'):
		
		return 280, 360, 0, 70

	else:
		
		print('Range does not exist.')

def Get_Mask_Range(Range):

	# Get model properties
	Model_Properties = Get_Model_Properties()
	Lon = Model_Properties['Lon']
	Lat = Model_Properties['Lat']

	Lon_Min, Lon_Max, Lat_Min, Lat_Max = Get_Range(Range)

	# Create mask array
	Mask_Range = (Lon>=Lon_Min)[None, :] & (Lon<=Lon_Max)[None, :] & (Lat>=Lat_Min)[:, None] & (Lat<=Lat_Max)[:, None]

	return Mask_Range

def Get_Mask_Month(Month, ENSO_Shifting=False):

	import numpy as np

	Mask_Month = np.full((12), False)
	
	if (isinstance(Month, str)):
	
		if (Month == 'All'):

			Mask_Month[:] = True
	
		elif (Month == 'SON'):

			Mask_Month[[8, 9, 10]] = True
	
		elif (Month == 'DJF'):

			Mask_Month[[11, 0, 1]] = True
	
		elif (Month == 'JA'):

			Mask_Month[[6, 7]] = True

	elif (isinstance(Month, list)):
		
		Mask_Month[Month] = True

	else:

		print('Get_Mask_Month: Error')
		quit()
	
	if (ENSO_Shifting):

		Mask_Month = np.roll(Mask_Month, -5)
	
	return Mask_Month

def Crop_Range(Data, Range, Range_Original='Global_Analysis'):

	import numpy as np

	# Get model properties
	Model_Properties = Get_Model_Properties()
	Lon = Model_Properties['Lon']
	Lat = Model_Properties['Lat']

	if (tuple(Data.shape[-2:]) != tuple([len(i) for i in Crop_Lat_Lon(Lat, Lon, Range_Original)])):

		print('Error in Crop_Range: given array does not meet original range.')
		return
	
	else:
		
		Lon_Min, Lon_Max, Lat_Min, Lat_Max = Get_Range(Range)
		Arg_Lon_Min = np.argmin(np.abs(Lon-Lon_Min))
		Arg_Lon_Max = np.argmin(np.abs(Lon-Lon_Max)) + 1
		Arg_Lat_Min = np.argmin(np.abs(Lat-Lat_Min))
		Arg_Lat_Max = np.argmin(np.abs(Lat-Lat_Max)) + 1
	
		if (Range_Original != 'Global_Analysis'):

			Lon_Min_Original, Lon_Max_Original, Lat_Min_Original, Lat_Max_Original = Get_Range(Range_Original)
			Arg_Lon_Min_Original = np.argmin(np.abs(Lon-Lon_Min_Original))
			Arg_Lon_Max_Original = np.argmin(np.abs(Lon-Lon_Max_Original)) + 1
			Arg_Lat_Min_Original = np.argmin(np.abs(Lat-Lat_Min_Original))
			Arg_Lat_Max_Original = np.argmin(np.abs(Lat-Lat_Max_Original)) + 1

			Arg_Lon_Min          = Arg_Lon_Min - Arg_Lon_Min_Original
			Arg_Lon_Max          = Arg_Lon_Max - Arg_Lon_Max_Original
			Arg_Lat_Min          = Arg_Lat_Min - Arg_Lat_Min_Original
			Arg_Lat_Max          = Arg_Lat_Max - Arg_Lat_Max_Original
		
		if ((Arg_Lat_Max in Lat[[0, -1]]) and (Arg_Lon_Max in Lon[[0, -1]])):

			Data = Data[..., Arg_Lat_Min:, Arg_Lon_Min:]
		
		elif (Arg_Lat_Max in Lat[[0, -1]]) and ~(Arg_Lon_Max in Lon[[0, -1]]):

			Data = Data[..., Arg_Lat_Min:, Arg_Lon_Min:Arg_Lon_Max]
		
		elif ~(Arg_Lat_Max in Lat[[0, -1]]) and (Arg_Lon_Max in Lon[[0, -1]]):

			Data = Data[..., Arg_Lat_Min:Arg_Lat_Max, Arg_Lon_Min:]
		
		else:

			Data = Data[..., Arg_Lat_Min:Arg_Lat_Max, Arg_Lon_Min:Arg_Lon_Max]
		
		return Data

def Crop_Lat_Lon(Lat, Lon, Range):
	
	import numpy as np

	# Get model properties
	Model_Properties = Get_Model_Properties()
	Lon = Model_Properties['Lon']
	Lat = Model_Properties['Lat']

	Lon_Min, Lon_Max, Lat_Min, Lat_Max = Get_Range(Range)
	Arg_Lon_Min = np.argmin(np.abs(Lon-Lon_Min))
	Arg_Lon_Max = np.argmin(np.abs(Lon-Lon_Max))
	Arg_Lat_Min = np.argmin(np.abs(Lat-Lat_Min))
	Arg_Lat_Max = np.argmin(np.abs(Lat-Lat_Max))

	Lat = Lat[Arg_Lat_Min:Arg_Lat_Max+1]
	Lon = Lon[Arg_Lon_Min:Arg_Lon_Max+1]

	return Lat, Lon

def Get_Mask_LandFraction(Type):

	import numpy as np

	# Get model properties
	Model_Properties = Get_Model_Properties()
	LandFraction = np.squeeze(Model_Properties['LandFraction'])
	
	if (Type == 'Land'):

		Mask_LandFraction = np.where(LandFraction>=0.5, LandFraction, 0)

	elif (Type == 'Ocean'):

		Mask_LandFraction = np.where(LandFraction<=0.5, 1-LandFraction, 0)
	
	elif (Type == 'Ocean_All'):

		Mask_LandFraction = np.where(LandFraction==0.0, 1, 0)
	
	elif (Type == 'All'):

		Mask_LandFraction = np.ones_like(LandFraction)

	return Mask_LandFraction

def Get_Weights_Lat(Range='Global_Analysis'):

	import numpy as np

	# Get model properties
	Model_Properties = Get_Model_Properties()
	Lon = Model_Properties['Lon']
	Lat = Model_Properties['Lat']

	Weights_Lat = np.tile(np.cos(np.deg2rad(Lat))[:, None], (1, Lon.shape[0]))

	return Weights_Lat

def Get_ENSO_List(Run, ONI_Threshold=1.5):

	import pandas as pd

	ENSO_List = pd.read_csv('Output_Data/File_ENSO/File_ENSO.csv')
	ENSO_List = ENSO_List[ENSO_List['Run']==Run]['ENSO_Threshold{}'.format(str(ONI_Threshold).replace('.', ''))].to_numpy()
	
	return ENSO_List

def Get_Mask_ENSO(Run, ENSO, ONI_Threshold=1.5):
	
	import numpy as np

	if (ENSO == 'All'):

		Mask_ENSO = np.full_like(Get_ENSO_List(Run), True)

	else:

		Mask_ENSO = np.where(Get_ENSO_List(Run, ONI_Threshold)==ENSO, True, False)
	
	return Mask_ENSO

def Calc_SpatialAverage(Data, Range, **kwargs):

	LandFraction_Type = kwargs.get('LandFraction_Type', 'Land')
	Mask_Range        = kwargs.get('Mask_Range', True)

	import numpy as np
	import warnings

	warnings.simplefilter('ignore')

	if (Mask_Range):

		Weight_Range        = np.where(Get_Mask_Range(Range), 1, 0)
		Weight_LandFraction = np.where(Get_Mask_LandFraction(LandFraction_Type), 1, 0)
		Weight_Lat          = Get_Weights_Lat(Range)
		Weights             = Weight_Range * Weight_LandFraction * Weight_Lat

		Data_NanMasked      = np.ma.masked_array(Data, np.isnan(Data))
		Data_SpatialAverage = np.ma.average(Data_NanMasked, axis=(-2, -1), weights=np.broadcast_to(Weights, Data_NanMasked.shape))

	else:

		Weight_LandFraction = np.where(Get_Mask_LandFraction(LandFraction_Type), 1, 0)
		Weight_Lat          = Get_Weights_Lat(Range)
		Weights             = Crop_Range(Weight_LandFraction * Weight_Lat, Range)

		Data_NanMasked      = np.ma.masked_array(Data, np.isnan(Data))
		Data_SpatialAverage = np.ma.average(Data_NanMasked, axis=(-2, -1), weights=np.broadcast_to(Weights, Data_NanMasked.shape))

	if (isinstance(Data_SpatialAverage, np.float64)):

		return Data_SpatialAverage
	
	else:

		return Data_SpatialAverage.filled(np.nan)

def Calc_AnnualMean(Data, t_Axis=0):

	import numpy as np
	import warnings
	import skimage.measure

	warnings.simplefilter('ignore', category=RuntimeWarning)

	# Create new dimension tuple
	NewDim         = [1] * len(Data.shape)
	NewDim[t_Axis] = 12

	Data_AnnualMean = skimage.measure.block_reduce(Data, block_size=tuple(NewDim), func=np.nanmean, cval=np.nan)

	return Data_AnnualMean

def Calc_SeasonalCycle(Data, **kwargs):

	import numpy as np
	import warnings

	Running_Filter = kwargs.get('Running_Filter', True)
	Filter_Range   = kwargs.get('Filter_Range', [-40, 40])

	warnings.simplefilter('ignore', category=RuntimeWarning)

	if (Running_Filter):

		if (np.ndim(Data) != 1):

			print('if Running_Filter is turn-on, data must be one-dimensional')
			return

		n_Year = Data.shape[0] // 12
		SeasonalCycle = np.full((n_Year, 12), np.nan)

		for i_Year in np.arange(n_Year):
			
			Filter_Start = i_Year + Filter_Range[0]
			Filter_End   = i_Year + Filter_Range[1]

			if (Filter_Start < 0): Filter_Start = 0
			if (Filter_End > n_Year): Filter_End = n_Year

			SeasonalCycle[i_Year, :] = np.nanmean(Data[Filter_Start*12:Filter_End*12].reshape((-1, 12)), axis=0)
		
		SeasonalCycle = SeasonalCycle.ravel()
	
	else:

		if (np.ndim(Data) == 1):

			SeasonalCycle = np.nanmean(Data.reshape((Data.shape[0]//12, 12)), 0)
			SeasonalCycle = np.tile(SeasonalCycle, (Data[0]//12))

		else:
			
			SeasonalCycle = np.nanmean(Data.reshape(Data.shape[0]//12, 12, *Data.shape[1:]), 0)
			SeasonalCycle = np.tile(SeasonalCycle, (Data.shape[0]//12, *[1] * (np.ndim(Data)-1)))

	return SeasonalCycle

def Calc_CI(Data, CI_Axis=0, Bool_return_Std=False):

	import numpy as np

	if (isinstance(Data, list)):

		# Confident between two samples with different size
		Mean_1 = np.nanmean(Data[0], axis=CI_Axis)
		Mean_2 = np.nanmean(Data[1], axis=CI_Axis)
		s_1 = np.nanstd(Data[0], axis=CI_Axis)
		s_2 = np.nanstd(Data[1], axis=CI_Axis)
		n_1 = np.sum(~np.isnan(Data[0]), axis=CI_Axis)
		n_2 = np.sum(~np.isnan(Data[1]), axis=CI_Axis)

		Sp  = np.sqrt(((n_1-1)*s_1**2+(n_2-1)*s_2**2)/(n_1+n_2-2))

		# Calculate the mean of seasonal cycle
		Data_Mean = Mean_1 - Mean_2

		# Calculate the CI of seasonal cycle
		Data_n    = np.min(np.array([n_1, n_2]))
		Data_CI   = [\
			Data_Mean - 1.96 * Sp * np.sqrt(1/n_1+1/n_2), \
			Data_Mean + 1.96 * Sp * np.sqrt(1/n_1+1/n_2), \
		]

		Data_Std = Sp

	else:
		
		# Calculate the mean of seasonal cycle
		Data_Mean = np.nanmean(Data, axis=CI_Axis)

		# Calculate the CI of seasonal cycle
		Data_Std  = np.nanstd(Data, axis=CI_Axis)
		Data_n    = np.sum(~np.isnan(Data), axis=CI_Axis)
		Data_CI   = [\
			Data_Mean - 1.96 * Data_Std / np.sqrt(Data_n), \
			Data_Mean + 1.96 * Data_Std / np.sqrt(Data_n), \
		]
	
	if (Bool_return_Std):
	
		return Data_Mean, Data_CI, Data_Std, Data_n

	else:

		return Data_Mean, Data_CI, Data_n