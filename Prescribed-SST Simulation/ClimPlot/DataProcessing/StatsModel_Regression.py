# StatsModel_Regression.py

# Import Module
import numpy as np
import pandas as pd

# Get models
def Get_RegressionModels(**kwargs):
	ModelType  = kwargs.get('ModelType', 'LinearRegression')
	n_jobs     = kwargs.get('n_jobs', 3)

	RegressionModels_Name = []
	RegressionModels_Object = []

	# Append linear regression
	import sklearn.linear_model
	RegressionModels_Name.append('LinearRegression')
	RegressionModels_Object.append(sklearn.linear_model.LinearRegression(n_jobs=n_jobs))
	
	if (ModelType!='LinearRegression'):

		# Append support vector regression
		import sklearn.svm
		RegressionModels_Name.append('SVR')
		RegressionModels_Object.append(sklearn.svm.SVR(kernel='poly', degree=2))

	return pd.DataFrame({'Name': RegressionModels_Name, 'Object': RegressionModels_Object})
