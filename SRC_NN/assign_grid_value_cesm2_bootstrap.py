from pandas import Series as se
from pandas import DataFrame as df
from scipy.io import loadmat
import scipy.stats
import pandas as pd
import numpy as np

import datetime

import tensorflow as tf
from tensorflow.keras import backend as K
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, Conv2D, Dropout, Flatten, MaxPooling2D
from tensorflow.keras.layers import Dropout
from tensorflow.keras.optimizers import SGD

from matplotlib import pyplot as plt

import random


########################################################
# Environment Setting
########################################################
is_server = 1
model_name = 'cesm2_clm5_mic_vr_v22' # 'cesm2_clm5_cen_vr_sce' # 'cesm2_clm5_cen_vr_v2'
time_domain = 'whole_time' # 'whole_time' 'before_1985' 'after_1985' 'random_half_1' 'random_half_2'

is_median_scaled = 0

for bootstrap_num in np.arange(1, 201):
	print('processing boottstrap ' + str(bootstrap_num))
	nn_training_name = 'exp_pc_cesm2_24' + '_bootstrap_' + str(bootstrap_num)
	# nn_training_name = 'exp_pc_cesm2_24' + '_cross_valid_0_' + str(bootstrap_num)
	########################################################
	# Data Path
	########################################################
	if is_server == 1:
		data_dir_output = '/GFPS8p/cess11/taof/datahub/ensemble/output_data/'
		data_dir_input = '/GFPS8p/cess11/taof/datahub/ensemble/input_data/'
	else:
		data_dir_output = '/Users/phoenix/Google_Drive/Tsinghua_Luo/Projects/DATAHUB/ENSEMBLE/OUTPUT_DATA/'
		data_dir_input = '/Users/phoenix/Google_Drive/Tsinghua_Luo/Projects/DATAHUB/ENSEMBLE/INPUT_DATA/'
	
	########################################################
	# Import basic info of environmental covariates 
	########################################################
	# variables used in training the NN
	# variables used in training the NN
	var4nn = ['Lon', 'Lat', \
	'ESA_Land_Cover', \
	# 'IGBP', \
	# 'Climate', \
	# 'Soil_Type', \
	# 'NPPmean', 'NPPmax', 'NPPmin', \
	# 'Veg_Cover', \
	'BIO1', 'BIO2', 'BIO3', 'BIO4', 'BIO5', 'BIO6', 'BIO7', 'BIO8', 'BIO9', 'BIO10', 'BIO11', 'BIO12', 'BIO13', 'BIO14', 'BIO15', 'BIO16', 'BIO17', 'BIO18', 'BIO19', \
	'Abs_Depth_to_Bedrock', \
	'Bulk_Density_0cm', 'Bulk_Density_30cm', 'Bulk_Density_100cm',\
	'CEC_0cm', 'CEC_30cm', 'CEC_100cm', \
	'Clay_Content_0cm', 'Clay_Content_30cm', 'Clay_Content_100cm', \
	'Coarse_Fragments_v_0cm', 'Coarse_Fragments_v_30cm', 'Coarse_Fragments_v_100cm', \
	# 'Depth_Bedrock_R', \
	'Garde_Acid', \
	'Occurrence_R_Horizon', \
	'pH_Water_0cm', 'pH_Water_30cm', 'pH_Water_100cm', \
	'Sand_Content_0cm', 'Sand_Content_30cm', 'Sand_Content_100cm', \
	'Silt_Content_0cm', 'Silt_Content_30cm', 'Silt_Content_100cm', \
	'SWC_v_Wilting_Point_0cm', 'SWC_v_Wilting_Point_30cm', 'SWC_v_Wilting_Point_100cm', \
	'Texture_USDA_0cm', 'Texture_USDA_30cm', 'Texture_USDA_100cm', \
	'USDA_Suborder', \
	'WRB_Subgroup', \
	# 'Drought', \
	'Elevation', \
	# 'Max_Depth'
	'Koppen_Climate_2018', \
	'cesm2_npp', 'cesm2_npp_std', \
	# 'cesm2_gpp', 'cesm2_gpp_std', \
	'cesm2_vegc', \
	'nbedrock']
	
	########################################################
	# Load trained NN model
	########################################################
	# define a joint loss
	para_mse = 1000*0.5
	para_ratio = 50*0.5
	def joint_loss (y_true, y_pred):
		# mse
		mse_loss = K.mean(K.square(y_true - y_pred))
		# mean absolute ratio  error
		ratio_loss = K.mean(K.abs((y_true - y_pred)/y_true))
		# return the joint loss
		return para_mse*mse_loss * para_ratio*ratio_loss
	
	def ratio_loss (y_true, y_pred):
		# mean absolute ratio  error
		ratio_loss = K.mean(K.abs((y_true - y_pred)/y_true))
		return ratio_loss
	
	model = tf.keras.models.load_model(data_dir_output + 'neural_networking/trained_model_' + model_name + '_' + time_domain + '_' + nn_training_name + '.h5', custom_objects={'joint_loss': joint_loss})
	
	########################################################
	# Load Global Pxiels
	########################################################
	# environmental info of global grids 
	grid_env_info_names = [\
	'Lon', 'Lat', 'Date', \
	'Rmean', 'Rmax', 'Rmin', \
	'ESA_Land_Cover', \
	'ET', \
	'IGBP', 'Climate', 'Soil_Type', 'NPPmean', 'NPPmax', 'NPPmin', \
	'Veg_Cover', \
	'BIO1', 'BIO2', 'BIO3', 'BIO4', 'BIO5', 'BIO6', 'BIO7', 'BIO8', 'BIO9', 'BIO10', 'BIO11', 'BIO12', 'BIO13', 'BIO14', 'BIO15', 'BIO16', 'BIO17', 'BIO18', 'BIO19', \
	'Abs_Depth_to_Bedrock', \
	'Bulk_Density_0cm', 'Bulk_Density_30cm', 'Bulk_Density_100cm',\
	'CEC_0cm', 'CEC_30cm', 'CEC_100cm', \
	'Clay_Content_0cm', 'Clay_Content_30cm', 'Clay_Content_100cm', \
	'Coarse_Fragments_v_0cm', 'Coarse_Fragments_v_30cm', 'Coarse_Fragments_v_100cm', \
	'Depth_Bedrock_R', \
	'Garde_Acid', \
	'Occurrence_R_Horizon', \
	'pH_Water_0cm', 'pH_Water_30cm', 'pH_Water_100cm', \
	'Sand_Content_0cm', 'Sand_Content_30cm', 'Sand_Content_100cm', \
	'Silt_Content_0cm', 'Silt_Content_30cm', 'Silt_Content_100cm', \
	'SWC_v_Wilting_Point_0cm', 'SWC_v_Wilting_Point_30cm', 'SWC_v_Wilting_Point_100cm', \
	'Texture_USDA_0cm', 'Texture_USDA_30cm', 'Texture_USDA_100cm', \
	'USDA_Suborder', \
	'WRB_Subgroup', \
	'Drought', \
	'Elevation', \
	'Max_Depth', \
	'Koppen_Climate_2018', \
	'cesm2_npp', 'cesm2_npp_std', \
	'cesm2_gpp', 'cesm2_gpp_std', \
	'cesm2_vegc', \
	'nbedrock']
	
	if is_median_scaled == 1:
		grid_env_info = loadmat(data_dir_input + 'data4nn/world_grid_envinfo_present_cesm2_clm5_cen_vr_v2_' + time_domain + '_median_scaled.mat')
	else:
		grid_env_info = loadmat(data_dir_input + 'data4nn/world_grid_envinfo_present_cesm2_clm5_cen_vr_v2_' + time_domain + '_maxmin_scaled.mat')
	
	
	grid_env_info = df(grid_env_info['grid_env_info'])
	
	grid_env_info.columns = grid_env_info_names
	
	# eliminate nan value    
	current_grid = np.array(grid_env_info.loc[:, var4nn])
	valid_grid_loc = np.array(range(len(current_grid[:, 1]))) + 1
	valid_grid_loc = np.reshape(valid_grid_loc, (len(current_grid[:, 1]), 1))
	
	for ivar in range(len(var4nn)):
		valid_grid_loc = valid_grid_loc[(np.isnan(current_grid[:, ivar]) == False)]
		current_grid = current_grid[(np.isnan(current_grid[:, ivar]) == False), :]
	
	## predict parameter values of each pixcel
	grid_predict = model.predict(current_grid)
		
	np.savetxt(data_dir_output + 'neural_networking/grid_para_result_' + model_name + '_' + time_domain + '_' + nn_training_name + '.csv', grid_predict, delimiter = ',')
	np.savetxt(data_dir_output + 'neural_networking/valid_grid_loc_' + model_name + '_' + time_domain + '_' + nn_training_name + '.csv', valid_grid_loc, delimiter = ',')
	
#end for bootstrap_num



