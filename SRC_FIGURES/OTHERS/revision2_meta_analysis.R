
library(R.matlab)
library(lme4)
library(lmerTest)

rm(list = ls())

setwd('/Users/phoenix/Google_Drive/R')

##################################################
## load emperical data
##################################################
grid_var_names = c('Lon', 'Lat', 'Date', 
                   'Rmean', 'Rmax', 'Rmin', 
                   'ESA_Land_Cover', 
                   'ET',
                   'IGBP', 'Climate', 'Soil_Type', 'NPPmean', 'NPPmax', 'NPPmin',
                   'Veg_Cover', 
                   'Annual Mean Temperature', 'Mean Diurnal Range', 'Isothermality', 'Temperature Seasonality', 'Max Temperature of Warmest Month', 'Min Temperature of Coldest Month', 'Temperature Annual Range', 'Mean Temperature of Wettest Quarter', 'Mean Temperature of Driest Quarter', 'Mean Temperature of Warmest Quarter', 'Mean Temperature of Coldest Quarter', 'Annual Precipitation', 'Precipitation of Wettest Month', 'Precipitation of Driest Month', 'Precipitation Seasonality', 'Precipitation of Wettest Quarter', 'Precipitation of Driest Quarter', 'Precipitation of Warmest Quarter', 'Precipitation of Coldest Quarter', 
                   'Abs_Depth_to_Bedrock',
                   'Bulk_Density_0cm', 'Bulk_Density_30cm', 'Bulk_Density_100cm',
                   'CEC_0cm', 'CEC_30cm', 'CEC_100cm',
                   'Clay_Content_0cm', 'Clay_Content_30cm', 'Clay_Content_100cm',
                   'Coarse_Fragments_v_0cm', 'Coarse_Fragments_v_30cm', 'Coarse_Fragments_v_100cm', 
                   'Depth_Bedrock_R', 
                   'Garde_Acid', 
                   'Occurrence_R_Horizon', 
                   'pH_Water_0cm', 'pH_Water_30cm', 'pH_Water_100cm', 
                   'Sand_Content_0cm', 'Sand_Content_30cm', 'Sand_Content_100cm',
                   'Silt_Content_0cm', 'Silt_Content_30cm', 'Silt_Content_100cm', 
                   'SWC_v_Wilting_Point_0cm', 'SWC_v_Wilting_Point_30cm', 'SWC_v_Wilting_Point_100cm', 
                   'Texture_USDA_0cm', 'Texture_USDA_30cm', 'Texture_USDA_100cm', 
                   'USDA_Suborder', 
                   'WRB_Subgroup', 
                   'Drought',
                   'Elevation',
                   'Max_Depth', 
                   'Koppen_Climate_2018', 
                   'cesm2_npp', 'cesm2_npp_std',
                   'cesm2_gpp', 'cesm2_gpp_std',
                   'cesm2_vegc',
                   'nbedrock')

var_nn_list = c('Bulk_Density_0cm', 'Bulk_Density_30cm', 'Bulk_Density_100cm',
                'Clay_Content_0cm', 'Clay_Content_30cm', 'Clay_Content_100cm',
                'Coarse_Fragments_v_0cm', 'Coarse_Fragments_v_30cm', 'Coarse_Fragments_v_100cm',
                'Sand_Content_0cm', 'Sand_Content_30cm', 'Sand_Content_100cm',
                'Silt_Content_0cm', 'Silt_Content_30cm', 'Silt_Content_100cm',
                'SWC_v_Wilting_Point_0cm', 'SWC_v_Wilting_Point_30cm', 'SWC_v_Wilting_Point_100cm', 
                'CEC_0cm', 'CEC_30cm', 'CEC_100cm',
                'Garde_Acid',
                'pH_Water_0cm', 'pH_Water_30cm', 'pH_Water_100cm',
                'Texture_USDA_0cm', 'Texture_USDA_30cm', 'Texture_USDA_100cm',
                'USDA_Suborder',
                'WRB_Subgroup',
                'Annual Mean Temperature', 'Mean Diurnal Range', 'Isothermality', 'Temperature Seasonality', 'Max Temperature of Warmest Month', 'Min Temperature of Coldest Month', 'Temperature Annual Range', 'Mean Temperature of Wettest Quarter', 'Mean Temperature of Driest Quarter', 'Mean Temperature of Warmest Quarter', 'Mean Temperature of Coldest Quarter', 'Annual Precipitation', 'Precipitation of Wettest Month', 'Precipitation of Driest Month', 'Precipitation Seasonality', 'Precipitation of Wettest Quarter', 'Precipitation of Driest Quarter', 'Precipitation of Warmest Quarter', 'Precipitation of Coldest Quarter', 'Koppen_Climate_2018',
                'ESA_Land_Cover', 
                'cesm2_npp', 'cesm2_npp_std', 
                'cesm2_vegc',
                'Lon', 'Lat',
                'Abs_Depth_to_Bedrock',
                'Occurrence_R_Horizon',
                'nbedrock',
                'Elevation')



soil_var_texture_list = c('Bulk_Density_0cm', 'Bulk_Density_30cm', 'Bulk_Density_100cm',
                          'Clay_Content_0cm', 'Clay_Content_30cm', 'Clay_Content_100cm',
                          'Coarse_Fragments_v_0cm', 'Coarse_Fragments_v_30cm', 'Coarse_Fragments_v_100cm',
                          'Sand_Content_0cm', 'Sand_Content_30cm', 'Sand_Content_100cm',
                          'Silt_Content_0cm', 'Silt_Content_30cm', 'Silt_Content_100cm',
                          'SWC_v_Wilting_Point_0cm', 'SWC_v_Wilting_Point_30cm', 'SWC_v_Wilting_Point_100cm')


soil_var_chemical_list = c('CEC_0cm', 'CEC_30cm', 'CEC_100cm',
                           'Garde_Acid',
                           'pH_Water_0cm', 'pH_Water_30cm', 'pH_Water_100cm',
                           'Texture_USDA_0cm', 'Texture_USDA_30cm', 'Texture_USDA_100cm',
                           'USDA_Suborder',
                           'WRB_Subgroup')

climate_var_list = c('Annual Mean Temperature', 'Mean Diurnal Range', 'Isothermality', 'Temperature Seasonality', 'Max Temperature of Warmest Month', 'Min Temperature of Coldest Month', 'Temperature Annual Range', 'Mean Temperature of Wettest Quarter', 'Mean Temperature of Driest Quarter', 'Mean Temperature of Warmest Quarter', 'Mean Temperature of Coldest Quarter', 'Annual Precipitation', 'Precipitation of Wettest Month', 'Precipitation of Driest Month', 'Precipitation Seasonality', 'Precipitation of Wettest Quarter', 'Precipitation of Driest Quarter', 'Precipitation of Warmest Quarter', 'Precipitation of Coldest Quarter', 'Koppen_Climate_2018')

vegetation_var_list =  c('ESA_Land_Cover', 
                         'cesm2_npp', 'cesm2_npp_std', 
                         'cesm2_vegc')

geography_var_list =  c('Lon', 'Lat',
                        'Abs_Depth_to_Bedrock',
                        'Occurrence_R_Horizon',
                        'nbedrock',
                        'Elevation')


env_info = readMat('/Users/phoenix/Google_Drive/Tsinghua_Luo/Projects/DATAHUB/CUE_Synthesis/cue_meta_envinfo.mat')
env_info = env_info$EnvInfo[ , 4:80]
colnames(env_info) = grid_var_names
cue_meta_info = readMat('/Users/phoenix/Google_Drive/Tsinghua_Luo/Projects/DATAHUB/CUE_Synthesis/cue_meta.mat')
cue_meta_info = cue_meta_info$meta.site.info
colnames(cue_meta_info) = c('site_id', 'lon', 'lat', 'mat', 'incubation_temp', 'map', 'depth', 'cue', 'mic', 'soc', 'pH', 'cue_method', 'source_id')

cue_meta_info[is.na(cue_meta_info[, 'mat']) == 1, 'mat'] = env_info[is.na(cue_meta_info[, 'mat']) == 1, 'Annual Mean Temperature']
cue_meta_info[is.na(cue_meta_info[, 'map']) == 1, 'map'] = env_info[is.na(cue_meta_info[, 'map']) == 1, 'Annual Precipitation']
cue_meta_info[is.na(cue_meta_info[, 'pH']) == 1, 'pH'] = env_info[is.na(cue_meta_info[, 'pH']) == 1, 'pH_Water_0cm']/10
cue_meta_info[cue_meta_info[ , 'cue_method'] < 3, 'cue_method'] = 1

valid_loc = which(cue_meta_info[ , 'source_id'] < 1000)
# valid_loc = which(cue_meta_info[ , 'source_id'] > 0)

current_data = data.frame(cue_meta_info[valid_loc, c('site_id', 'source_id', 'mat', 'map', 'incubation_temp', 'cue', 'mic', 'soc', 'depth', 'pH', 'cue_method')])
current_data$non_mic = current_data$soc - current_data$mic/1000

cor.test(current_data[ , c('cue')], current_data[ , c('mic')])
cor.test(current_data[ , c('cue')], current_data[ , c('soc')])
cor.test(current_data[ , c('cue')], current_data[ , c('non_mic')])
cor.test(current_data[ , c('mic')], current_data[ , c('non_mic')])

####################################################################
# CUE-SOC different method
####################################################################
current_data = data.frame(cue_meta_info[valid_loc, c('site_id', 'source_id', 'mat', 'cue', 'soc', 'depth', 'cue_method')])
current_data$site_id = as.factor(current_data$site_id)
current_data$source_id = as.factor(current_data$source_id)
current_data$cue_method = as.factor(current_data$cue_method)


mix_model_isotope_c = lmer(log10(soc) ~ cue + depth + mat + (1 | source_id), data = current_data[current_data$cue_method == 1, ])
mix_model_isotope_o = lmer(log10(soc) ~ cue + depth + mat + (1 | source_id), data = current_data[current_data$cue_method == 3, ])

r_squared_isotope_c = summary(lm(current_data$soc[current_data$cue_method == 1] ~ predict(mix_model_isotope_c)))
r_squared_isotope_c 

r_squared_isotope_o = summary(lm(current_data$soc[current_data$cue_method == 3] ~ predict(mix_model_isotope_o)))
r_squared_isotope_o 

summary(mix_model_isotope_c)
summary(mix_model_isotope_o)

####################################################################
# CUE-SOC (log value)
####################################################################
current_data = data.frame(cue_meta_info[valid_loc, c('site_id', 'source_id', 'mat', 'map', 'incubation_temp', 'cue', 'soc', 'depth', 'pH', 'cue_method')])
current_data$site_id = as.factor(current_data$site_id)
current_data$source_id = as.factor(current_data$source_id)
current_data$cue_method = as.factor(current_data$cue_method)

cor.test(current_data[ , c('cue')], current_data[ , c('depth')])
cor.test(current_data[ , c('cue')], current_data[ , c('mat')])
cor.test(current_data[ , c('mat')], current_data[ , c('depth')])

# 
# summary(lm(current_data$cue ~ current_data$mat + current_data$pH))
# 
# ggplot(data = current_data) +
#   geom_point(aes(x = cue, y = soc, group = cue_method, color = cue_method)) +
#   geom_smooth(aes(x = cue, y = soc, group = cue_method, color = cue_method), fill = NA,  method = 'lm') +
#   scale_y_continuous(trans = 'log10')
# 
## mixed effects model

mix_model_1 = lmer(log10(soc) ~ cue + depth + mat + (1 | source_id), data = current_data)
summary(mix_model_1)
car::vif(mix_model_1)

r_squared_mixed = summary(lm(log10(current_data$soc) ~ predict(mix_model_1)))
r_squared_mixed 

fix_effect_model = summary(mix_model_1)
fix_effect_pred = fix_effect_model$coefficients[1, 1] + 
  fix_effect_model$coefficients[2, 1]*current_data$cue + 
  fix_effect_model$coefficients[3, 1]*current_data$depth + 
  fix_effect_model$coefficients[4, 1]*current_data$mat
r_squared_fixed = summary(lm(log10(current_data$soc) ~ fix_effect_pred))
r_squared_fixed


#------------------
mix_model_2 = lmer(log10(soc) ~ cue + depth + mat + mat*cue + (1 | source_id), data = current_data)
summary(mix_model_2)
car::vif(mix_model_2)

r_squared_mixed = summary(lm(log10(current_data$soc) ~ predict(mix_model_2)))
r_squared_mixed 

fix_effect_model = summary(mix_model_2)
fix_effect_pred = fix_effect_model$coefficients[1, 1] + 
  fix_effect_model$coefficients[2, 1]*current_data$cue + 
  fix_effect_model$coefficients[3, 1]*current_data$depth + 
  fix_effect_model$coefficients[4, 1]*current_data$mat
r_squared_fixed = summary(lm(log10(current_data$soc) ~ fix_effect_pred))
r_squared_fixed

#------------------
mix_model_3 = lmer(log10(soc) ~ cue + depth + mat + mat*cue + (cue | source_id), data = current_data)
summary(mix_model_3)
car::vif(mix_model_3)

r_squared_mixed = summary(lm(current_data$soc ~ predict(mix_model_3)))
r_squared_mixed 

fix_effect_model = summary(mix_model_3)
fix_effect_pred = fix_effect_model$coefficients[1, 1] + 
  fix_effect_model$coefficients[2, 1]*current_data$cue + 
  fix_effect_model$coefficients[3, 1]*current_data$depth + 
  fix_effect_model$coefficients[4, 1]*current_data$mat
r_squared_fixed = summary(lm(log10(current_data$soc) ~ fix_effect_pred))
r_squared_fixed


anova(mix_model_1, mix_model_2, mix_model_3, refit = FALSE)

####################################################################
# CUE-SOC (original value)
####################################################################
current_data = data.frame(cue_meta_info[valid_loc, c('site_id', 'source_id', 'mat', 'map', 'incubation_temp', 'cue', 'soc', 'depth', 'pH', 'cue_method')])
current_data$site_id = as.factor(current_data$site_id)
current_data$source_id = as.factor(current_data$source_id)
current_data$cue_method = as.factor(current_data$cue_method)

# ggplot(data = current_data) +
#   geom_point(aes(x = cue, y = soc, group = cue_method, color = cue_method)) +
#   geom_smooth(aes(x = cue, y = soc, group = cue_method, color = cue_method), fill = NA,  method = 'lm') +
#   scale_y_continuous(trans = 'log10')
# 
## mixed effects model

mix_model_1 = lmer(soc ~ cue + depth + mat + (1 | source_id), data = current_data)
summary(mix_model_1)
car::vif(mix_model_1)

r_squared_mixed = summary(lm(current_data$soc ~ predict(mix_model_1)))
r_squared_mixed 

fix_effect_model = summary(mix_model_1)
fix_effect_pred = fix_effect_model$coefficients[1, 1] + 
  fix_effect_model$coefficients[2, 1]*current_data$cue + 
  fix_effect_model$coefficients[3, 1]*current_data$depth + 
  fix_effect_model$coefficients[4, 1]*current_data$mat
r_squared_fixed = summary(lm(current_data$soc ~ fix_effect_pred))
r_squared_fixed


mix_model_2 = lmer(soc ~ cue + depth + mat + (cue | source_id), data = current_data)
summary(mix_model_2)
car::vif(mix_model_2)

r_squared_mixed = summary(lm(current_data$soc ~ predict(mix_model_2)))
r_squared_mixed 

fix_effect_model = summary(mix_model_2)
fix_effect_pred = fix_effect_model$coefficients[1, 1] + 
  fix_effect_model$coefficients[2, 1]*current_data$cue + 
  fix_effect_model$coefficients[3, 1]*current_data$depth + 
  fix_effect_model$coefficients[4, 1]*current_data$mat
r_squared_fixed = summary(lm(current_data$soc ~ fix_effect_pred))
r_squared_fixed


mix_model_3 = lmer(soc ~ cue + depth + mat + cue*mat + (1 | source_id), data = current_data)
summary(mix_model_3)
car::vif(mix_model_3)

r_squared_mixed = summary(lm(current_data$soc ~ predict(mix_model_3)))
r_squared_mixed 

fix_effect_model = summary(mix_model_3)
fix_effect_pred = fix_effect_model$coefficients[1, 1] + 
  fix_effect_model$coefficients[2, 1]*current_data$cue + 
  fix_effect_model$coefficients[3, 1]*current_data$depth + 
  fix_effect_model$coefficients[4, 1]*current_data$mat
r_squared_fixed = summary(lm(current_data$soc ~ fix_effect_pred))
r_squared_fixed


anova(mix_model_1, mix_model_2, mix_model_3, refit = FALSE)
