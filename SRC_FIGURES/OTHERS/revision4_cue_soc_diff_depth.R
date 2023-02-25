## Packages
library(R.matlab)
library(maps)
library(ggplot2)
library(R.matlab)
library(cowplot)
library(jcolors)
library(viridis)
library(cowplot)
library(scales)

library(lme4)
library(lmerTest)
library(corrplot)

library(ppcor)
# dev.off()
##
rm(list = ls())

setwd('/Users/phoenix/Google_Drive/R')

#############################################################################
# microbial data path
#############################################################################
version_num = 22
model_name = paste('cesm2_clm5_mic_vr_v', version_num, sep = '')
model_name_beta = paste(strsplit(model_name, split = '_')[[1]], collapse = '.')

data_dir_output = '/Users/phoenix/Google_Drive/Tsinghua_Luo/Projects/DATAHUB/ENSEMBLE/OUTPUT_DATA/'
data_dir_input = '/Users/phoenix/Google_Drive/Tsinghua_Luo/Projects/DATAHUB/ENSEMBLE/INPUT_DATA/'

var_names = c('ProfileNum', 'ProfileID', 'LayerNum', 'Lon', 'Lat', 'Date',
              'Rmean', 'Rmax', 'Rmin',
              'ESA_Land_Cover',
              'ET', 
              'IGBP', 'Climate', 'Soil_Type', 'NPPmean', 'NPPmax', 'NPPmin',
              'Veg_Cover', 
              'BIO1', 'BIO2', 'BIO3', 'BIO4', 'BIO5', 'BIO6', 'BIO7', 'BIO8', 'BIO9', 'BIO10', 'BIO11', 'BIO12', 'BIO13', 'BIO14', 'BIO15', 'BIO16', 'BIO17', 'BIO18', 'BIO19', 
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
              'nbedrock',
              'R_Squared')

#----------------------------------------------------
# microbial profile info
#----------------------------------------------------
profile_env_info = readMat(paste(data_dir_input, '/wosis_2019_snap_shot/wosis_2019_snapshot_hugelius_mishra_env_info.mat', sep = ''))
profile_env_info = profile_env_info$EnvInfo
colnames(profile_env_info) = var_names

# load climate info
global_climate = profile_env_info[ , 'Koppen_Climate_2018']
var_climate = array(NA, dim = c(length(global_climate), 1))
# koppen climatte class 2
var_climate = global_climate
var_climate[which(global_climate == 1)] = 101 # Af
var_climate[which(global_climate == 2)] = 102 # Am
var_climate[which(global_climate == 3)] = 103 # Aw
var_climate[which(global_climate >= 4 & global_climate <= 5)] = 104 # BwX
var_climate[which(global_climate >= 6 & global_climate <= 7)] = 105 # BsW
var_climate[which(global_climate >= 8 & global_climate <= 10)] = 106 # CsX
var_climate[which(global_climate >= 11 & global_climate <= 13)] = 107 # CwX
var_climate[which(global_climate >= 14 & global_climate <= 16)] = 108 # CfX
var_climate[which(global_climate >= 17 & global_climate <= 20)] = 109 # DsX
var_climate[which(global_climate >= 21 & global_climate <= 24)] = 110 # DwX
var_climate[which(global_climate >= 25 & global_climate <= 28)] = 111 # DfX
var_climate[which(global_climate >= 29 & global_climate <= 30)] = 112 # E


#----------------------------------------------------
# microbial load da summary
#----------------------------------------------------
# parameters
para_name = c('bio', 'cryo', 'q10', 'w_scaling', 
              'tau4cwd', 'tau4l1', 'tau4l2', 'tau4s1', 'tau4s2_death', 'tau4s2_enz', 'tau4s3', 'tau4s4', 
              'mm_const_assim', 'mm_const_decom', 
              'fcwdl2', 
              'pl1s1', 'pl2s1', 'pl3s4', 'l1_cue', 'l2l3_cue', 
              'mic_cue', 'pdeath2soc',
              'beta')
npara = length(para_name)
da_summary_mic_para = readMat(paste(data_dir_output, '/mcmc_summary_', model_name, '/', model_name, '_para_mean.mat', sep = ''))
da_summary_mic_para = da_summary_mic_para$para.mean

colnames(da_summary_mic_para) = para_name

# data assimilation summary
col_list_relation_summary = c('profile_num', 'mic_cue', 
                              'soc_0_30cm', 'soc_30_100cm', 'soc_100_200cm', 'mic_0_30cm', 
                              'mic_stock', 'enz_stock', 'doc_stock', 'poc_stock', 'soc_total', 
                              'bulk_A', 'bulk_I', 'bulk_K', 'bulk_V', 'bulk_Xi', 'bulk_E', 'bulk_NPP')
relation_summary = readMat(paste(data_dir_output, '/mcmc_summary_', model_name, '/', model_name, '_da_summary_mic.mat', sep = ''))
relation_summary = relation_summary$da.summary.mic
colnames(relation_summary) = col_list_relation_summary

mm_ratio_assim = (da_summary_mic_para[ , 'mm_const_assim']*(3000 - 300) + 300)/relation_summary[ , 'doc_stock']/0.3
mm_ratio_decom = (da_summary_mic_para[ , 'mm_const_decom']*(10**6 - 10**5) + 10**5)/relation_summary[ , 'soc_0_30cm']/0.3

relation_summary[ , 'soc_0_30cm'] = relation_summary[ , 'soc_0_30cm']/0.3/profile_env_info[ , 'Bulk_Density_0cm'] # unit gc/kg
relation_summary[ , 'soc_30_100cm'] = relation_summary[ , 'soc_30_100cm']/0.7/profile_env_info[ , 'Bulk_Density_30cm'] # unit gc/kg
relation_summary[ , 'soc_100_200cm'] = relation_summary[ , 'soc_100_200cm']/1/profile_env_info[ , 'Bulk_Density_100cm'] # unit gc/kg

relation_summary[ , 'mic_0_30cm'] = relation_summary[ , 'mic_0_30cm']/0.3/profile_env_info[ , 'Bulk_Density_0cm'] # unit gc/kg

#----------------------------------------------------
# microbial load bulk CUE
#----------------------------------------------------
mic_bulk_process_summary = relation_summary[ , c('bulk_A', 'bulk_I', 'bulk_K', 'bulk_V', 'bulk_Xi', 'bulk_E', 'bulk_NPP')]

valid_profile_loc = which(is.na(relation_summary[ , 'mic_cue']) == 0)

#############################################################################
# microbial figure
#############################################################################
#----------------------------------------------------
# michaelis menten ratios histograms
#----------------------------------------------------
current_data = data.frame(rbind(cbind(mm_ratio_assim[valid_profile_loc], 1), 
                                cbind(mm_ratio_decom[valid_profile_loc], 2)))
colnames(current_data) = c('ratio', 'class')

color_scheme = c('#DC3220', '#005AB5')

jpeg(paste('./Ensemble/revision4_michaelis_menten_ratio_', model_name, '.jpeg', sep = ''), width = 10, height = 10, units = 'in', res = 300)

ggplot() +
  geom_histogram(data = current_data, aes(x = ratio, y = ..density.., color = as.factor(class), fill = as.factor(class)), position = 'identity', alpha = 0.3, na.rm = TRUE, size = 1, bins = 30) + 
  geom_density(data = current_data, aes(x = ratio, color = as.factor(class)), fill = NA, alpha = 0.5, na.rm = TRUE, size = 3) + 
  geom_vline(xintercept = c(100), linetype = 'longdash', size = 2, color = 'grey') + 
  scale_x_continuous(trans = 'log10', n.breaks = 6, labels = trans_format('log10', math_format(10^.x))) +
  scale_fill_manual(name = '', values = color_scheme, labels = c('Assimilation', 'Decomposition')) +
  scale_color_manual(name = '', values = color_scheme, labels = c('Assimilation', 'Decomposition')) +
  # change the background to black and white
  theme_classic() + 
  # add title
  labs(title = '', x = 'Km/[Substrate]', y = 'Density') + 
  # change the legend properties
  theme(legend.text = element_text(size = 30), legend.title = element_text(size = 30))  +
  theme(legend.justification = c(1, 1), legend.position = c(1, 1.05), legend.background = element_rect(fill = NA)) + 
  theme(legend.key = element_rect(color = NA, fill = NA), legend.key.size = unit(0.5, 'inch')) +
  # modify the position of title
  theme(plot.title = element_text(hjust = 0.5, size = 50)) + 
  # modify the font size
  # modify the margin
  # theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
  theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), 'inch')) +
  theme(axis.text=element_text(size = 30), axis.title = element_text(size = 35), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1), axis.ticks.length = unit(0.12, 'inch')) 

dev.off()

#----------------------------------------------------
# cue soc relation 0 - 30 cm
#----------------------------------------------------
current_data = data.frame(cbind(relation_summary[valid_profile_loc, c('soc_0_30cm')],
                                relation_summary[valid_profile_loc, c('mic_cue')],
                                valid_profile_loc, 
                                var_climate[valid_profile_loc],
                                profile_env_info[valid_profile_loc, c('BIO1')]))

colnames(current_data) = c('soc', 'cue', 'profile_num', 'climate', 'mat')

current_data = current_data[is.na(apply(current_data, 1, sum, na.rm = FALSE)) == 0, ]

#---------------mix model---------------#
mix_model = lmer(log10(soc) ~ cue + (cue | climate), data = current_data)
mix_model_summary_1 = summary(mix_model)
mix_model_summary_1
fit_function_1 = function(x) {10**(mix_model_summary_1$coefficients[1, 1] + x*(mix_model_summary_1$coefficients[2, 1]))}

r_squared_mixed = summary(lm(log10(current_data$soc[is.na(current_data$soc) == 0]) ~ predict(mix_model)))
r_squared_mixed 

if (mix_model_summary_1$coefficients[1, 5] < 0.001) {
  p_intercept = paste(' < 0.001', sep = '')
} else {
  p_intercept = paste(' = ', round(mix_model_summary_1$coefficients[1, 5], 3), sep = '')
}

if (mix_model_summary_1$coefficients[2, 5] < 0.001) {
  p_slope = paste(' < 0.001', sep = '')
} else {
  p_slope = paste(' = ', round(mix_model_summary_1$coefficients[2, 5], 3), sep = '')
}

text_data = data.frame('x_axis' = c(0.0), 
                       'y_axis' = c(1), 
                       'equation' = paste('log10(SOC) = ', round(mix_model_summary_1$coefficients[1, 1], 2), ' + ', round(mix_model_summary_1$coefficients[2, 1], 2), '*CUE', '\nP(Intercept)', p_intercept, ', P(CUE)', p_slope, '\nExplained Variation = ', round(r_squared_mixed$r.squared*100, 2), '%', sep = ''))


color_scheme = c('#b30000', '#e34a33', '#fc8d59', '#fdbb84', '#fdd49e', '#fef0d9', '#f0f9e8', '#ccebc5', '#a8ddb5', '#7bccc4', '#43a2ca', '#0868ac')


p_mic_cue_soc_0_30cm =
  ggplot(data = current_data) + 
  geom_point(data = current_data, aes(x = cue, y = soc, color = mat), size = 1, alpha = 0.7, na.rm = TRUE) +
  geom_function(fun = fit_function_1, size = 2, color = 'black') + 
  geom_text(data = text_data, aes(x = x_axis, y = y_axis, label = equation), vjust = 1, hjust = 0, size = 8.5, show.legend = FALSE) + 
  scale_color_gradientn(name = 'MAT (°C)', colors = rev(color_scheme), limits = c(-5, 30), oob = scales::squish) +
  scale_x_continuous(trans = 'identity') +
  scale_y_continuous(limits = c(0.1, 1000), trans = 'log10', labels = trans_format('log10', math_format(10^.x))) + 
  # change the background to black and white
  theme_classic() + 
  # add title
  labs(title = '', x = 'Carbon use efficiency', y = expression(paste('Soil organic carbon (g C kg'^'-1', ')', sep = ''))) + 
  # change the legend properties
  guides(color = guide_colorbar(direction = 'horizontal', barwidth = 15, barheight = 2.5, title.position = 'right', title.hjust = 0, title.vjust = 0.8, label.hjust = 0.2, frame.linewidth = 0), reverse = FALSE) +
  theme(legend.text = element_text(size = 25), legend.title = element_text(size = 25))  +
  theme(legend.justification = c(0, 1), legend.position = c(0, 1), legend.background = element_rect(fill = NA)) + 
  # modify the position of title
  theme(plot.title = element_text(hjust = 0.5, size = 50)) + 
  # modify the font size
  # modify the margin
  # theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
  theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), 'inch')) +
  theme(axis.text=element_text(size = 30), axis.title = element_text(size = 35), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1), axis.ticks.length = unit(0.12, 'inch')) 

#----------------------------------------------------
# cue soc relation 30 - 100 cm
#----------------------------------------------------
current_data = data.frame(cbind(relation_summary[valid_profile_loc, c('soc_30_100cm')],
                                relation_summary[valid_profile_loc, c('mic_cue')],
                                valid_profile_loc, 
                                var_climate[valid_profile_loc],
                                profile_env_info[valid_profile_loc, c('BIO1')]))

colnames(current_data) = c('soc', 'cue', 'profile_num', 'climate', 'mat')

current_data = current_data[is.na(apply(current_data, 1, sum, na.rm = FALSE)) == 0, ]
#---------------mix model---------------#
mix_model = lmer(log10(soc) ~ cue + (cue | climate), data = current_data)
mix_model_summary_2 = summary(mix_model)
mix_model_summary_2
fit_function_2 = function(x) {10**(mix_model_summary_2$coefficients[1, 1] + x*(mix_model_summary_2$coefficients[2, 1]))}

r_squared_mixed = summary(lm(log10(current_data$soc[is.na(current_data$soc) == 0]) ~ predict(mix_model)))
r_squared_mixed 


if (mix_model_summary_2$coefficients[1, 5] < 0.001) {
  p_intercept = paste(' < 0.001', sep = '')
} else {
  p_intercept = paste(' = ', round(mix_model_summary_2$coefficients[1, 5], 3), sep = '')
}

if (mix_model_summary_2$coefficients[2, 5] < 0.001) {
  p_slope = paste(' < 0.001', sep = '')
} else {
  p_slope = paste(' = ', round(mix_model_summary_2$coefficients[2, 5], 3), sep = '')
}

text_data = data.frame('x_axis' = c(0.0), 
                       'y_axis' = c(1000), 
                       'equation' = paste('log10(SOC) = ', round(mix_model_summary_2$coefficients[1, 1], 2), ' + ', round(mix_model_summary_2$coefficients[2, 1], 2), '*CUE', '\nP(Intercept)', p_intercept, ', P(CUE)', p_slope, '\nExplained Variation = ', round(r_squared_mixed$r.squared*100, 2), '%', sep = ''))


p_mic_cue_soc_30_100cm =
  ggplot(data = current_data) + 
  geom_point(data = current_data, aes(x = cue, y = soc, color = mat), size = 1, alpha = 0.7, na.rm = TRUE) +
  geom_function(fun = fit_function_2, size = 2, color = 'black') + 
  # geom_smooth(data = current_data, aes(x = cue, y = soc_0_30cm), method = 'lm', color = 'black', fill = 'black', size = 2) +
  geom_text(data = text_data, aes(x = x_axis, y = y_axis, label = equation), vjust = 1, hjust = 0, size = 8.5, show.legend = FALSE) + 
  scale_color_gradientn(name = 'MAT (°C)', colors = rev(color_scheme), limits = c(-5, 30), oob = scales::squish) +
  scale_x_continuous(trans = 'identity') +
  scale_y_continuous(limits = c(0.1, 1000), trans = 'log10', labels = trans_format('log10', math_format(10^.x))) + 
  # change the background to black and white
  theme_classic() + 
  # add title
  labs(title = '', x = 'Carbon use efficiency', y = expression(paste('Soil organic carbon (g C kg'^'-1', ')', sep = ''))) + 
  # change the legend properties
  guides(color = guide_colorbar(direction = 'horizontal', barwidth = 15, barheight = 2.5, title.position = 'top', title.hjust = 0, label.hjust = 0, frame.linewidth = 0), reverse = FALSE) +
  theme(legend.text = element_text(size = 25), legend.title = element_text(size = 25))  +
  theme(legend.justification = c(0, 0), legend.position = 'None', legend.background = element_rect(fill = NA)) + 
  # modify the position of title
  theme(plot.title = element_text(hjust = 0.5, size = 50)) + 
  # modify the font size
  # modify the margin
  # theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
  theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), 'inch')) +
  theme(axis.text=element_text(size = 30), axis.title = element_text(size = 35), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1), axis.ticks.length = unit(0.12, 'inch')) 

#----------------------------------------------------
# cue soc relation 100 - max cm
#----------------------------------------------------
current_data = data.frame(cbind(relation_summary[valid_profile_loc, c('soc_100_200cm')],
                                relation_summary[valid_profile_loc, c('mic_cue')],
                                valid_profile_loc, 
                                var_climate[valid_profile_loc],
                                profile_env_info[valid_profile_loc, c('BIO1')]))

colnames(current_data) = c('soc', 'cue', 'profile_num', 'climate', 'mat')


current_data = current_data[is.na(apply(current_data, 1, sum, na.rm = FALSE)) == 0, ]

#---------------mix model---------------#
mix_model = lmer(log10(soc) ~ cue + (cue | climate), data = current_data)
mix_model_summary_3 = summary(mix_model)
mix_model_summary_3
fit_function_3 = function(x) {10**(mix_model_summary_3$coefficients[1, 1] + x*(mix_model_summary_3$coefficients[2, 1]))}

r_squared_mixed = summary(lm(log10(current_data$soc[is.na(current_data$soc) == 0]) ~ predict(mix_model)))
r_squared_mixed 


if (mix_model_summary_3$coefficients[1, 5] < 0.001) {
  p_intercept = paste(' < 0.001', sep = '')
} else {
  p_intercept = paste(' = ', round(mix_model_summary_3$coefficients[1, 5], 3), sep = '')
}

if (mix_model_summary_3$coefficients[2, 5] < 0.001) {
  p_slope = paste(' < 0.001', sep = '')
} else {
  p_slope = paste(' = ', round(mix_model_summary_3$coefficients[2, 5], 3), sep = '')
}

text_data = data.frame('x_axis' = c(0.0), 
                       'y_axis' = c(1000), 
                       'equation' = paste('log10(SOC) = ', round(mix_model_summary_3$coefficients[1, 1], 2), ' + ', round(mix_model_summary_3$coefficients[2, 1], 2), '*CUE', '\nP(Intercept)', p_intercept, ', P(CUE)', p_slope, '\nExplained Variation = ', round(r_squared_mixed$r.squared*100, 2), '%', sep = ''))


p_mic_cue_soc_100_200cm =
  ggplot(data = current_data) + 
  geom_point(data = current_data, aes(x = cue, y = soc, color = mat), size = 1, alpha = 0.7, na.rm = TRUE) +
  geom_function(fun = fit_function_3, size = 2, color = 'black') + 
  # geom_smooth(data = current_data, aes(x = cue, y = soc_0_30cm), method = 'lm', color = 'black', fill = 'black', size = 2) +
  geom_text(data = text_data, aes(x = x_axis, y = y_axis, label = equation), vjust = 1, hjust = 0, size = 8.5, show.legend = FALSE) + 
  scale_color_gradientn(name = 'MAT (°C)', colors = rev(color_scheme), limits = c(-5, 30), oob = scales::squish) +
  scale_x_continuous(trans = 'identity') +
  scale_y_continuous(limits = c(0.01, 1000), trans = 'log10', labels = trans_format('log10', math_format(10^.x))) + 
  # change the background to black and white
  theme_classic() + 
  # add title
  labs(title = '', x = 'Carbon use efficiency', y = expression(paste('Soil organic carbon (g C kg'^'-1', ')', sep = ''))) + 
  # change the legend properties
  guides(color = guide_colorbar(direction = 'horizontal', barwidth = 15, barheight = 2.5, title.position = 'top', title.hjust = 0, label.hjust = 0, frame.linewidth = 0), reverse = FALSE) +
  theme(legend.text = element_text(size = 25), legend.title = element_text(size = 25))  +
  theme(legend.justification = c(0, 0), legend.position = 'None', legend.background = element_rect(fill = NA)) + 
  # modify the position of title
  theme(plot.title = element_text(hjust = 0.5, size = 50)) + 
  # modify the font size
  # modify the margin
  # theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
  theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), 'inch')) +
  theme(axis.text=element_text(size = 30), axis.title = element_text(size = 35), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1), axis.ticks.length = unit(0.12, 'inch')) 
  
########################################
# final figure
########################################
jpeg(paste('./Ensemble/revision4_cue_soc_relation_depth.jpeg', sep = ''), width = 24, height = 8, units = 'in', res = 300)
p_mic_cue_soc_0_30cm = p_mic_cue_soc_0_30cm + 
  labs(title = '0 - 30cm')
p_mic_cue_soc_30_100cm = p_mic_cue_soc_30_100cm + 
  labs(title = '30 - 100cm', y = '  ')
p_mic_cue_soc_100_200cm  = p_mic_cue_soc_100_200cm + 
  labs(title = '> 100cm',y = '  ')

plot_grid(  p_mic_cue_soc_0_30cm, p_mic_cue_soc_30_100cm, p_mic_cue_soc_100_200cm,
            labels = c('a', 'b', 'c'),
            label_size = 50,
            label_x = 0.07, label_y = 1.02,
            label_fontfamily = 'Arial',
            label_fontface = 'bold',
            nrow = 1)
dev.off()

