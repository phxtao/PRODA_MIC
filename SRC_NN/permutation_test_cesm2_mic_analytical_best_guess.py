from pandas import Series as se
from pandas import DataFrame as df
from scipy.io import loadmat
import scipy.stats
import pandas as pd
import numpy as np

from netCDF4 import Dataset

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
is_server = 0

model_name = 'cesm2_clm5_mic_vr_v22'

time_domain = 'whole_time' # 'whole_time', 'before_1985', 'after_1985'

nn_training_name = 'exp_pc_cesm2_23'

cross_valid_num = 1
# nn_training_name = 'exp_pc_cesm2_23' + '_cross_valid_' + str(cross_valid_num)

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
# Import para info after bayesian method
########################################################
# laod para info after MCMC
para_without_trans = loadmat(data_dir_output + 'mcmc_summary_' + model_name + '/' + model_name + '_para_mean.mat')
para_without_trans = df(para_without_trans['para_mean'])

para_names = ['bio', 'cryo', 'q10', 'w_scaling', 'tau4cwd', 'tau4l1', 'tau4l2', 'tau4s1', 'tau4s2_death', 'tau4s2_enz', 'tau4s3', 'tau4s4', 'mm_const_assim', 'mm_const_decom', 'fcwdl2', 'pl1s1', 'pl2s1', 'pl3s4', 'l1_cue', 'l2l3_cue', 'mic_cue', 'pdeath2soc', 'beta']

para_names_cate = ['V', 'V', 'Xi', 'Xi', 'K', 'K', 'K', 'A', 'K', 'K', 'K', 'K', 'A', 'K', 'E', 'E', 'E', 'E', 'A', 'A', 'A', 'E', 'I']

# index order A K Xi V I E
para_names_index = [4, 4, 3, 3, 2, 2, 2, 1, 2, 2, 2, 2, 1, 2, 6, 6, 6, 6, 1, 1, 1, 6, 5]

para_without_trans.columns = para_names

# environmental info of soil profiles
env_info = loadmat(data_dir_input + 'wosis_2019_snap_shot/wosis_2019_snapshot_hugelius_mishra_env_info.mat')

env_info = env_info['EnvInfo']

col_max_min = loadmat(data_dir_input + 'data4nn/world_grid_envinfo_present_cesm2_clm5_cen_vr_v2_whole_time_col_max_min.mat')
col_max_min = col_max_min['col_max_min']
for ivar in np.arange(3, len(col_max_min[:, 0])):
	env_info[:, ivar] = (env_info[:, ivar] - col_max_min[ivar, 0])/(col_max_min[ivar, 1] - col_max_min[ivar, 0])
	env_info[(env_info[:, ivar] > 1), ivar] = 1
	env_info[(env_info[:, ivar] < 0), ivar] = 0

env_info = df(env_info)

env_info_names = ['ProfileNum', 'ProfileID', 'LayerNum', 'Lon', 'Lat', 'Date', \
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
'nbedrock', \
'R_Squared']

env_info.columns = env_info_names

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
# 'Max_Depth', \
'Koppen_Climate_2018', \
'cesm2_npp', 'cesm2_npp_std', \
# 'cesm2_gpp', 'cesm2_gpp_std', \
'cesm2_vegc', \
'nbedrock', \
]

soil_var_structure =  ['Clay_Content_0cm', 'Clay_Content_30cm', 'Clay_Content_100cm', \
'Sand_Content_0cm', 'Sand_Content_30cm', 'Sand_Content_100cm', \
'Silt_Content_0cm', 'Silt_Content_30cm', 'Silt_Content_100cm', \
'Texture_USDA_0cm', 'Texture_USDA_30cm', 'Texture_USDA_100cm', \
'USDA_Suborder', \
'WRB_Subgroup', \
'Bulk_Density_0cm', 'Bulk_Density_30cm', 'Bulk_Density_100cm',\
'Coarse_Fragments_v_0cm', 'Coarse_Fragments_v_30cm', 'Coarse_Fragments_v_100cm', \
'SWC_v_Wilting_Point_0cm', 'SWC_v_Wilting_Point_30cm', 'SWC_v_Wilting_Point_100cm']


soil_var_chemical = ['CEC_0cm', 'CEC_30cm', 'CEC_100cm', \
'Garde_Acid', \
'pH_Water_0cm', 'pH_Water_30cm', 'pH_Water_100cm']

climate_var = ['BIO1', 'BIO2', 'BIO3', 'BIO4', 'BIO5', 'BIO6', 'BIO7', 'BIO8', 'BIO9', 'BIO10', 'BIO11', 'BIO12', 'BIO13', 'BIO14', 'BIO15', 'BIO16', 'BIO17', 'BIO18', 'BIO19', \
'Koppen_Climate_2018'\
]

vegetation_var =  ['ESA_Land_Cover', \
'cesm2_npp', 'cesm2_npp_std', \
'cesm2_vegc' \
]

geography_var =  ['Lon', 'Lat', \
'Abs_Depth_to_Bedrock', \
'Occurrence_R_Horizon', \
'nbedrock', \
'Elevation' \
]

########################################################
# Load Training Information
########################################################
category_num = 5
corr_importance_cate_matric = np.empty((category_num, 2*len(para_names)))
corr_importance_cate_matric[:, :] = None

mse_importance_cate_matric = np.empty((category_num, 2*len(para_names)))
mse_importance_cate_matric[:, :] = None


corr_importance_matric = np.empty((len(var4nn), 2*len(para_names)))
corr_importance_matric[:, :] = None

mse_importance_matric = np.empty((len(var4nn), 2*len(para_names)))
mse_importance_matric[:, :] = None


permu_size = 1000

########################################################
# Permutation by category to category 
########################################################
para_cate_num = len(np.unique(para_names_cate))
corr_summary = np.empty((permu_size*cross_valid_num, para_cate_num))
corr_summary[:, :] = None
mse_summary = np.empty((permu_size*cross_valid_num, para_cate_num))
mse_summary[:, :] = None

for ivar in range(category_num):	
	for icross_valid in range(cross_valid_num):
		nn_training_name = 'exp_pc_cesm2_23' + '_cross_valid_0_' + str(7)
		
		nn_site_loc = pd.read_csv(data_dir_output + 'neural_networking/nn_site_loc_' + model_name + '_' + time_domain + '_' + nn_training_name + '.csv', sep  = ',', header = None)
		nn_site_loc = (nn_site_loc).loc[:, 0]
		
		
		para_without_trans.index = env_info.loc[:, 'ProfileNum']
		env_info.index = env_info.loc[:, 'ProfileNum']
		
		test_y = np.array(para_without_trans.loc[nn_site_loc, :])
		test_x = np.array(env_info.loc[nn_site_loc, var4nn])
		
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
		
		# Standard test as benchmark
		nn_predict_benckmark = model.predict(test_x)
	
		corr_para_benckmark = [None]*para_cate_num
		mse_para_benckmark = [None]*para_cate_num
		for ipara_cate in range(para_cate_num):
			obs_middle = np.concatenate(test_y[:, para_names_index == np.unique(para_names_index)[ipara_cate]])
			mod_middle = np.concatenate(nn_predict_benckmark[:, para_names_index == np.unique(para_names_index)[ipara_cate]])
			corr_para_benckmark[ipara_cate] = np.corrcoef(obs_middle, mod_middle)[0, 1]
			mse_para_benckmark[ipara_cate] = sum((obs_middle - mod_middle)**2)/len(obs_middle)
	
		# permutation
		for ipermu in range(permu_size):
			print('Processing permutation on category ' + str(ivar)  + ' cross valid ' + str(icross_valid) + ' permulation time ' + str(ipermu))
			# assign random value to a certain categorical env. var.
			test_x_permu = np.empty((len(test_y), len(var4nn)))
			test_x_permu[:, :] = None
			test_x_permu[:, :] = test_x
			test_x_permu = df(test_x_permu)
			test_x_permu.columns = var4nn
			if ivar == 0:
				random_envinfo = np.random.rand(1)*np.ones((len(test_y), len(soil_var_structure)))
				test_x_permu.loc[:, soil_var_structure] = test_x_permu.loc[:, soil_var_structure]*0 + random_envinfo*1
				test_x_permu = np.asarray(test_x_permu)
				test_x_permu[test_x_permu < 0] = 0
				test_x_permu[test_x_permu > 1] = 1
			if ivar == 1:
				random_envinfo = np.random.rand(1)*np.ones((len(test_y), len(soil_var_chemical)))
				test_x_permu.loc[:, soil_var_chemical] = test_x_permu.loc[:, soil_var_chemical]*0 + random_envinfo*1
				test_x_permu = np.asarray(test_x_permu)
				test_x_permu[test_x_permu < 0] = 0
				test_x_permu[test_x_permu > 1] = 1
			if ivar == 2:
				random_envinfo = np.random.rand(1)*np.ones((len(test_y), len(climate_var)))
				test_x_permu.loc[:, climate_var] =  test_x_permu.loc[:, climate_var]*0 + random_envinfo*1
				test_x_permu = np.asarray(test_x_permu)
				test_x_permu[test_x_permu < 0] = 0
				test_x_permu[test_x_permu > 1] = 1
			if ivar == 3:
				random_envinfo = np.random.rand(1)*np.ones((len(test_y), len(vegetation_var)))
				test_x_permu.loc[:, vegetation_var] = test_x_permu.loc[:, vegetation_var]*0 + random_envinfo*1
				test_x_permu = np.asarray(test_x_permu)
				test_x_permu[test_x_permu < 0] = 0
				test_x_permu[test_x_permu > 1] = 1
			if ivar == 4:
				random_envinfo = np.random.rand(1)*np.ones((len(test_y), len(geography_var)))
				test_x_permu.loc[:, geography_var] = test_x_permu.loc[:, geography_var]*0 + random_envinfo*1
				test_x_permu = np.asarray(test_x_permu)
				test_x_permu[test_x_permu < 0] = 0
				test_x_permu[test_x_permu > 1] = 1
				
			# do prediction after permutation
			nn_predict_permu = model.predict(test_x_permu)
			# calculate the importance of permutated var.
			for ipara_cate in range(para_cate_num):
				obs_middle = np.concatenate(test_y[:, para_names_index == np.unique(para_names_index)[ipara_cate]])
				mod_middle = np.concatenate(nn_predict_permu[:, para_names_index == np.unique(para_names_index)[ipara_cate]])
				corr_summary[ipermu + icross_valid*permu_size, ipara_cate] = np.corrcoef(obs_middle, mod_middle)[0, 1] / corr_para_benckmark[ipara_cate]
				mse_summary[ipermu + icross_valid*permu_size, ipara_cate] = np.nansum((obs_middle - mod_middle)**2)/len(obs_middle) / mse_para_benckmark[ipara_cate]
	# save summary	
	np.savetxt(data_dir_output + 'neural_networking/permutation_test_mse_summary_cate2cate_best_guess_var_' + str(ivar) + '_' + model_name + '_' + time_domain + '_' + nn_training_name + '.csv', mse_summary, delimiter = ',')
	
	# summarize the importance
	for ipara_cate in range(para_cate_num):
		mean_corr_error = np.nanmean(corr_summary[:, ipara_cate])
		std_corr_error = np.nanstd(corr_summary[:, ipara_cate])
		mean_mse_error = np.nanmean(mse_summary[:, ipara_cate])
		std_mse_error = np.nanstd(mse_summary[:, ipara_cate])
		
		corr_importance_cate_matric[ivar, ipara_cate] = mean_corr_error
		corr_importance_cate_matric[ivar, ipara_cate+para_cate_num] = std_corr_error
		mse_importance_cate_matric[ivar, ipara_cate] = mean_mse_error
		mse_importance_cate_matric[ivar, ipara_cate+para_cate_num] = std_mse_error
	
# np.savetxt(data_dir_output + 'neural_networking/permutation_test_corr_cate2cate_best_guess_' + model_name + '_' + time_domain + '_' + nn_training_name + '.csv', corr_importance_cate_matric, delimiter = ',')
np.savetxt(data_dir_output + 'neural_networking/permutation_test_mse_cate2cate_best_guess_' + model_name + '_' + time_domain + '_' + nn_training_name + '.csv', mse_importance_cate_matric, delimiter = ',')



