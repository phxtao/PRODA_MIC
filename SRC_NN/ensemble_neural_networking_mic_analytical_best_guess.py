import sys
from pandas import Series as se
from pandas import DataFrame as df
from scipy.io import loadmat
import scipy.stats
import pandas as pd
import numpy as np

import tensorflow as tf
from tensorflow.keras import backend as K
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, Conv2D, Dropout, Flatten, MaxPooling2D
from tensorflow.keras.layers import Dropout
from tensorflow.keras.layers import Activation
from tensorflow.keras.optimizers import SGD
from tensorflow.keras.utils import get_custom_objects
from tensorflow.keras.regularizers import l2
from tensorflow.keras.callbacks import EarlyStopping
from tensorflow.keras.callbacks import ModelCheckpoint
from tensorflow.keras.models import load_model

from matplotlib import pyplot as plt

import random

########################################################
# Environment Setting
########################################################
seed = int(sys.argv[1])

random.seed(seed)

is_server = 1
model_name = 'cesm2_clm5_mic_vr_v22'
time_domain = 'whole_time' # 'whole_time', 'before_1985', 'after_1985', 'random_half_1', 'random_half_2'
cross_valid_num = 0

nn_training_name = 'exp_pc_cesm2_24' + '_cross_valid_' + str(cross_valid_num) + '_' + str(seed)

is_median_scaled = 0
nn_loss = 'joint_loss' # 'mean_squared_error'
nn_optimizer = 'adadelta'
nn_batch_size = 32 # 16  when sample size is small
nn_epochs = 1200*5 # 1200*2
early_stop_patience = 1200
nn_layer_num = [256, 512, 512, 256]
nn_drop_ratio = [0.2]*len(nn_layer_num) #[0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
nn_l2_regularizer = [0.0]*len(nn_layer_num)

nn_split_ratio = 0.1

use_custom_activation = 0
nn_activation = [None]*len(nn_layer_num)

if use_custom_activation == 1:
	# define activation function
	def custom_activation(x):
		custom_activation = tf.keras.activations.relu(x, alpha = 0.1)
		return custom_activation
	
	for ilayer in range(len(nn_layer_num)):
		get_custom_objects().update({'custom_activation_'+str(ilayer): Activation(custom_activation)})
		nn_activation[ilayer] = 'custom_activation_'+str(ilayer)
else:
	nn_activation = ['relu', 'relu', 'relu', 'relu']

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

stat_r2 = loadmat(data_dir_output + 'mcmc_summary_' + model_name + '/' + model_name + '_stat_r2.mat')
stat_r2 = df(stat_r2['stat_r2'])
stat_r2 = np.max(stat_r2, axis = 1)

stat_gr = loadmat(data_dir_output + 'mcmc_summary_' + model_name + '/' + model_name + '_para_gr.mat')
stat_gr = df(stat_gr['para_gr'])
stat_gr = np.mean(stat_gr, axis = 1)

para_names = ['bio', 'cryo', 'q10', 'w_scaling', \
'tau4cwd', 'tau4l1', 'tau4l2', 'tau4s1', 'tau4s2_death', 'tau4s2_enz', 'tau4s3', 'tau4s4', \
'mm_const_assim', 'mm_const_decom', \
'fcwdl2', \
'pl1s1', 'pl2s1', 'pl3s4', 'l1_cue', 'l2l3_cue', \
'mic_cue', 'pdeath2soc', \
'beta']

# para_names = ['bio', 'cryo', 'q10', 'w_scaling', \
# 'tau4cwd', 'tau4l1', 'tau4l2', 'tau4s1', 'tau4s2_death', 'tau4s2_enz', 'tau4s3', 'tau4s4', \
# 'mm_const_assim', 'mm_const_decom', \
# 'fcwdl2', \
# 'pl1s1', 'pl2s1', 'pl3s4', 'l1_cue', \
# 'mic_cue', 'pdeath2soc', \
# 'beta']

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
'nbedrock']

########################################################
# specify the profiles 
########################################################
clean_loc= np.where((stat_r2 > 0.0) & (stat_gr < 1.05))[0]

eligible_loc = loadmat(data_dir_input + 'data4nn/eligible_profile_loc_0_cesm2_clm5_cen_vr_v2_whole_time.mat')
eligible_loc = eligible_loc['eligible_loc_0'][:, 0] - 1

eligible_loc = np.intersect1d(clean_loc, eligible_loc)

valid_loc = env_info.loc[eligible_loc, 'ProfileNum']
currentdata_y = np.array(para_without_trans.loc[eligible_loc, :])
currentdata_x = np.array(env_info.loc[eligible_loc, var4nn])

valid_loc = valid_loc[(np.isnan(currentdata_y[:, 0]) == False)] 
currentdata_x = currentdata_x[(np.isnan(currentdata_y[:, 0]) == False), :]
currentdata_y = currentdata_y[(np.isnan(currentdata_y[:, 0]) == False), :]
# eliminate nan value
for ivar in range(len(var4nn)):
	currentdata_y = currentdata_y[(np.isnan(currentdata_x[:, ivar]) == False), :]
	valid_loc = valid_loc[(np.isnan(currentdata_x[:, ivar]) == False)]
	currentdata_x = currentdata_x[(np.isnan(currentdata_x[:, ivar]) == False), :]

train_loc = np.random.choice(np.arange(0, len(currentdata_x[:, 0])), size = round((1-nn_split_ratio)*len(currentdata_x[:, 0])), replace = False)
test_loc = np.setdiff1d(np.arange(0, len(currentdata_x[:, 0])), train_loc)
# train_loc = np.arange(0, len(currentdata_x[:, 0]))

########################################################
# Configuration NN
########################################################

# split into input and outputs
train_x = currentdata_x[train_loc, :]
train_y = currentdata_y[train_loc, :]
 
test_x = currentdata_x[test_loc, :]
test_y = currentdata_y[test_loc, :]

print(train_x.shape, train_y.shape)

# define a joint loss
para_mse = 1000*0.5
para_ratio = 50*0.5
def joint_loss (y_true, y_pred):
	# mse
	mse_loss = K.mean(K.square(y_true - y_pred))
	# mean absolute ratio  error
	ratio_loss = K.mean(K.abs((y_true - y_pred)/y_true))
	# return the joint loss
	return para_mse*mse_loss # * para_ratio*ratio_loss

def ratio_loss (y_true, y_pred):
	# mean absolute ratio  error
	ratio_loss = K.mean(K.abs((y_true - y_pred)/y_true))
	return ratio_loss


# design network

model = Sequential()

for ilayer in range(len(nn_layer_num)):
	if use_custom_activation == 1:
		if ilayer == 0:
			model.add(Dense(nn_layer_num[ilayer], input_dim = len(var4nn)))
			model.add(Activation(custom_activation, name = nn_activation[ilayer]))
			model.add(Dropout(nn_drop_ratio[ilayer]))
		else:
			model.add(Dense(nn_layer_num[ilayer]))
			model.add(Activation(custom_activation, name = nn_activation[ilayer]))
			model.add(Dropout(nn_drop_ratio[ilayer]))
	else:
		if ilayer == 0:
			model.add(Dense(nn_layer_num[ilayer], kernel_regularizer=l2(nn_l2_regularizer[ilayer]), input_dim = len(var4nn), activation=nn_activation[ilayer]))
			model.add(Dropout(nn_drop_ratio[ilayer]))
		else:
			model.add(Dense(nn_layer_num[ilayer], kernel_regularizer=l2(nn_l2_regularizer[ilayer]), activation=nn_activation[ilayer]))
			model.add(Dropout(nn_drop_ratio[ilayer]))

model.add(Dense(len(para_names), activation='linear'))

if nn_loss == 'joint_loss':
	model.compile(loss = joint_loss, optimizer = nn_optimizer, metrics = ['accuracy'])
else:
	model.compile(loss = nn_loss, optimizer = nn_optimizer, metrics = ['accuracy'])

model.summary()

# early stopping
early_stop = EarlyStopping(monitor = 'val_loss', mode = 'min', verbose = 2, patience = early_stop_patience)

model_check = ModelCheckpoint(data_dir_output + 'neural_networking/trained_model_' + model_name + '_' + time_domain + '_' + nn_training_name + '.h5', monitor = 'val_loss', mode = 'min', verbose = 2, save_best_only = True)

########################################################
# NN Operation
########################################################
# fit network
history = model.fit(x = train_x, y = train_y, epochs = nn_epochs, batch_size = nn_batch_size, validation_split = nn_split_ratio, verbose = 2, callbacks=[early_stop, model_check])
# load best model
best_model = load_model(data_dir_output + 'neural_networking/trained_model_' + model_name + '_' + time_domain + '_' + nn_training_name + '.h5', custom_objects={'joint_loss': joint_loss})

# fig = plt.figure()    

nn_predict = best_model.predict(train_x)
corr_para = [None]*len(para_names)
for ipara in range(len(para_names)):
	corr_para[ipara] = np.corrcoef(train_y[:, ipara], nn_predict[:, ipara])[0, 1]
print(corr_para)

test_nn_predict = best_model.predict(test_x)
######################################################
# Output
######################################################
nn_site_loc =  np.asarray(valid_loc)[test_loc]
np.savetxt(data_dir_output + 'neural_networking/nn_para_result_' + model_name + '_' + time_domain + '_' + nn_training_name + '.csv', test_nn_predict, delimiter = ',')
np.savetxt(data_dir_output + 'neural_networking/nn_site_loc_' + model_name + '_' + time_domain + '_' + nn_training_name + '.csv', nn_site_loc, delimiter = ',')

