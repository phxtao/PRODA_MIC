## Packages
library(R.matlab)
library(maps)
library(ggplot2)
library(R.matlab)
library(cowplot)
library(jcolors)
library(viridis)
library(cowplot)
library(lme4)
library(lmerTest)
library(scales)

# dev.off()
##
rm(list = ls())

setwd('/Users/phoenix/Google_Drive/R')

## Jet colorbar function
jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
diff.colors <- colorRampPalette(c("#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "white", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B"))

#############################################################################
# Data Path
#############################################################################
data_dir_output = '/Users/phoenix/Google_Drive/Tsinghua_Luo/Projects/DATAHUB/ENSEMBLE/OUTPUT_DATA/'
data_dir_input = '/Users/phoenix/Google_Drive/Tsinghua_Luo/Projects/DATAHUB/ENSEMBLE/INPUT_DATA/'

#############################################################################
# function to increase vertical spacing between legend keys
#############################################################################
# @clauswilke
draw_key_polygon3 <- function(data, params, size) {
  lwd <- min(data$size, min(size) / 4)
  
  grid::rectGrob(
    width = grid::unit(0.6, "npc"),
    height = grid::unit(0.6, "npc"),
    gp = grid::gpar(
      col = data$colour,
      fill = alpha(data$fill, data$alpha),
      lty = data$linetype,
      lwd = lwd * .pt,
      linejoin = "mitre"
    ))
}

# register new key drawing function, 
# the effect is global & persistent throughout the R session
GeomBar$draw_key = draw_key_polygon3


#################################################################################
# meta-analysis
#################################################################################
## load emperical data
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
cue_meta_info[is.na(cue_meta_info[, 'pH']) == 1, 'pH'] = env_info[is.na(cue_meta_info[, 'pH']) == 1, 'pH_Water_0cm']/10
cue_meta_info[cue_meta_info[ , 'cue_method'] < 3, 'cue_method'] = 1


#-----------------------------------------------
# cue comparison between proda and meta-analysis
#-----------------------------------------------
proda_predict_meta = readMat(paste(data_dir_output, 'mcmc_summary_cesm2_clm5_mic_vr_v22/proda_prediction_on_meta.mat', sep = ''))
proda_predict_meta = proda_predict_meta$da.summary.mic

colnames(proda_predict_meta) = c('profile_num', 'mic_cue', 
               'soc_0_30cm', 'soc_30_100cm', 'soc_100_200cm', 'mic_0_30cm', 
               'mic_stock', 'enz_stock', 'doc_stock', 'poc_stock', 'soc_total', 
               'bulk_A', 'bulk_I', 'bulk_K', 'bulk_V', 'bulk_Xi', 'bulk_E', 'bulk_NPP')


valid_loc = which(cue_meta_info[ , 'source_id'] < 1000)
#------------------inter source cue
current_data_inter_cue = cbind(cue_meta_info[valid_loc, c('cue', 'soc', 'depth', 'mat')], 
                               proda_predict_meta[ , 'mic_cue'], 
                               proda_predict_meta[ , 'soc_0_30cm']/env_info[valid_loc, 'Bulk_Density_0cm'], 
                               cue_meta_info[valid_loc, 'source_id'])
current_data_inter_cue = data.frame(current_data_inter_cue)
colnames(current_data_inter_cue) = c('meta_cue', 'meta_soc', 'depth', 'mat', 'proda_cue', 'proda_soc', 'source_id')


# mixed effects models meta_cue vs proda_cue
mix_model = lmer(meta_cue ~ proda_cue + (1 | source_id), data = current_data_inter_cue)
mix_model_summary_inter_cue = summary(mix_model)
mix_model_summary_inter_cue
r_squared_mixed = summary(lm(current_data_inter_cue$meta_cue ~ predict(mix_model)))
r_squared_mixed

# whetehr slope is different from 1
mix_model = lmer((meta_cue) ~ proda_cue + (1 | source_id), offset = proda_cue, data = current_data_inter_cue)
mix_model_summary_inter_cue = summary(mix_model)
mix_model_summary_inter_cue
anova(mix_model) 

cor.test(current_data_inter_cue$meta_cue, current_data_inter_cue$proda_cue)

summary(lm(current_data_inter_cue$meta_cue ~ current_data_inter_cue$proda_cue))



# mixed model 

current_data_inter_cue = current_data_inter_cue[is.na(apply(current_data_inter_cue, 1, sum, na.rm = FALSE)) == 0, ]

mix_model_meta = lmer(log10(meta_soc) ~ meta_cue + depth + mat + (1 | source_id), data = current_data_inter_cue)
mix_model_summary_meta = summary(mix_model_meta)
mix_model_summary_meta

mix_model_proda = lmer(log10(proda_soc) ~ proda_cue + mat + (1 | source_id), data = current_data_inter_cue)
mix_model_summary_proda = summary(mix_model_proda)
mix_model_summary_proda

anova(mix_model_proda, mix_model_meta)

# ggplot() + 
#   geom_abline(intercept = 0, slope = 1) + 
#   geom_point(data = current_data_inter_cue, aes(x = proda, y = meta)) + 
#   geom_smooth(data = current_data_inter_cue, aes(x = proda, y = meta), method = 'lm') +
#   # scale_x_continuous(limits = c(0, 0.8)) +
#   # scale_y_continuous(limits = c(0, 0.8)) +
#   theme_classic()


#---------------------------------------------------
# CUE-SOC relationship
#---------------------------------------------------

#---------------mix model microbial vs non-microbial biomass---------------#
current_data_meta = data.frame(cue_meta_info[valid_loc, c('site_id', 'source_id', 'mat', 'cue', 'mic', 'soc', 'depth', 'cue_method')])
current_data_meta = current_data_meta[is.na(apply(current_data_meta, 1, sum, na.rm = FALSE)) == 0, ]

mix_model = lmer((mic/1000) ~ cue + depth + mat + (1 | source_id), data = current_data_meta)
mix_model_summary_meta = summary(mix_model)
mix_model_summary_meta
r_squared_mixed = summary(lm(current_data_meta$mic/1000 ~ predict(mix_model)))
r_squared_mixed


mix_model = lmer((soc - mic/1000) ~ cue + depth + mat + (1 | source_id), data = current_data_meta)
mix_model_summary_meta = summary(mix_model)
mix_model_summary_meta
r_squared_mixed = summary(lm((current_data_meta$soc - current_data_meta$mic/1000) ~ predict(mix_model)))
r_squared_mixed

#---------------mix model---------------#
current_data_meta = data.frame(cue_meta_info[valid_loc, c('site_id', 'source_id', 'mat', 'cue', 'soc', 'depth', 'cue_method')])
current_data_meta = current_data_meta[is.na(apply(current_data_meta, 1, sum, na.rm = FALSE)) == 0, ]

mix_model = lmer(log10(soc) ~ cue + depth + mat + (1 | source_id), data = current_data_meta)
mix_model_summary_meta = summary(mix_model)
mix_model_summary_meta

fit_function_meta = function(x) {10**(mix_model_summary_meta$coefficients[1, 1] + x*(mix_model_summary_meta$coefficients[2, 1]))}

r_squared_mixed = summary(lm(log10(current_data_meta$soc) ~ predict(mix_model)))
r_squared_mixed

if (mix_model_summary_meta$coefficients[1, 5] < 0.001) {
  p_intercept = paste(' < 0.001', sep = '')
} else {
  p_intercept = paste(' = ', round(mix_model_summary_meta$coefficients[1, 5], 3), sep = '')
}

if (mix_model_summary_meta$coefficients[2, 5] < 0.001) {
  p_slope = paste(' < 0.001', sep = '')
} else {
  p_slope = paste(' = ', round(mix_model_summary_meta$coefficients[2, 5], 3), sep = '')
}

text_data = data.frame('x_axis' = c(0.01), 
                       'y_axis' = c(1000), 
                       'equation' = paste('log10(SOC) = ', round(mix_model_summary_meta$coefficients[1, 1], 2), ' + ', round(mix_model_summary_meta$coefficients[2, 1], 2), '*CUE', '\nP(Intercept)', p_intercept, ', P(CUE)', p_slope, '\nExplained Variation = ', round(r_squared_mixed$r.squared*100, 2), '%', sep = ''))

color_scheme = c('#005AB5', '#DC3220')

p_meta_cue_soc =
  ggplot() +
  geom_point(data = current_data_meta, aes(x = cue, y = soc, size = depth), color = 'black', alpha = 0.7, na.rm = TRUE) +
  geom_function(fun = fit_function_meta, size = 2, color = 'black') + 
  scale_x_continuous(trans = 'identity') + 
  scale_y_continuous(limits = c(0.1, 1000), trans = 'log10', labels = trans_format('log10', math_format(10^.x))) + 
  scale_size_continuous(name = 'Depth (cm)', range = c(3, 9), breaks = c(5, 15, 30)) + 
  # change the background to black and white
  geom_text(data = text_data, aes(x = x_axis, y = y_axis, label = equation), vjust = 1, hjust = 0, size = 9, show.legend = FALSE) + 
  # change the background to black and white
  theme_classic() +
  # add title
  labs(title = '', x = 'Carbon use efficiency', y = expression(paste('Soil organic carbon  (g C kg'^'-1', ')', sep = ''))) + 
  # modify the position of title
  theme(plot.title = element_text(hjust = 0.5, size = 50)) + 
  theme(legend.justification = c(0, 0), legend.position = c(0, 0), legend.background = element_rect(fill = NA), legend.text.align = 0) +
  theme(legend.text = element_text(size = 30), legend.title = element_text(size = 30))  +
  theme(legend.key = element_rect(color = NA, fill = NA), legend.key.size = unit(0.5, 'inch')) +
  theme(legend.direction = 'horizontal') + 
  # modify the margin
  # theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
  theme(plot.margin = unit(c(0., 0.2, 0.2, 0.2), 'inch')) +
  theme(axis.text=element_text(size = 30), axis.title = element_text(size = 35), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1), axis.ticks.length = unit(0.12, 'inch')) 

#---------------------------------------------------
# MAT-CUE relationship
#---------------------------------------------------
corr_test = cor.test(current_data_meta$mat, current_data_meta$cue)

if (corr_test$p.value < 0.001) { 
  text_data = data.frame('x_axis' = c(-3), 
                         'y_axis' = c(0.8), 
                         'equation' = paste('Pearson corelation = ', round(corr_test$estimate, 2), 
                                            '\nP-value < 0.001', 
                                            sep = ''))
} else {
  text_data = data.frame('x_axis' = c(-3), 
                         'y_axis' = c(0.8), 
                         'equation' = paste('Pearson corelation = ', round(corr_test$estimate, 2), 
                                            '\nP-value = ', round(corr_test$p.value, 3), 
                                            sep = ''))
}

p_meta_mat_cue =
  ggplot() +
  geom_point(data = current_data_meta, aes(x = mat, y = cue, size = depth), color = 'snow4', alpha = 0.7, na.rm = TRUE) +
  geom_smooth(data = current_data_meta, aes(x = mat, y = cue), size = 2, color = 'black', fill = 'black', method = 'lm', alpha = 0.2) +
  scale_x_continuous(trans = 'identity') +
  scale_y_continuous(limits = c(0, 0.8), trans = 'identity') +
  scale_size_continuous(name = 'Depth (cm)', range = c(3, 9), breaks = c(5, 15, 30)) +
  # change the background to black and white
  geom_text(data = text_data, aes(x = x_axis, y = y_axis, label = equation), vjust = 1, hjust = 0, size = 9, show.legend = FALSE) + 
  # change the background to black and white
  theme_classic() +
  # add title
  labs(title = '', y = 'Carbon use efficiency', x = 'Mean annual temperature (°C)') + 
  # modify the position of title
  theme(plot.title = element_text(hjust = 0.5, size = 50)) + 
  theme(legend.justification = c(0, 0), legend.position = c(0, 0), legend.background = element_rect(fill = NA), legend.text.align = 0) +
  theme(legend.text = element_text(size = 30), legend.title = element_text(size = 30))  +
  theme(legend.key = element_rect(color = NA, fill = NA), legend.key.size = unit(0.5, 'inch')) +
  theme(legend.direction = 'horizontal') + 
  # modify the margin
  # theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
  theme(plot.margin = unit(c(0., 0.2, 0.2, 0.2), 'inch')) +
  theme(axis.text=element_text(size = 30), axis.title = element_text(size = 35), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1), axis.ticks.length = unit(0.12, 'inch')) 



#---------------------------------------------------
# CUE by different method
#---------------------------------------------------
color_scheme = c('#DC3220', '#005AB5')

jpeg(paste('./Ensemble/revision4_cue_diff_method.jpeg', sep = ''), width = 10, height = 10, units = 'in', res = 300)

ggplot() +
  geom_histogram(data = current_data_meta, aes(x = cue, y = ..density.., color = as.factor(cue_method), fill = as.factor(cue_method)), position = 'identity', alpha = 0.3, na.rm = TRUE, size = 1, bins = 20) + 
  scale_x_continuous(n.breaks = 6) +
  scale_fill_manual(name = '', values = color_scheme, labels = c(expression(paste(''^'13', 'C/'^'14', 'C', sep = '')), expression(paste(''^'18', 'O')))) +
  scale_color_manual(name = '', values = color_scheme, c(expression(paste(''^'13', 'C/'^'14', 'C', sep = '')), expression(paste(''^'18', 'O')))) +
  # change the background to black and white
  theme_classic() + 
  # add title
  labs(title = '', x = 'Microbial carbon use efficiency', y = 'Density') + 
  # change the legend properties
  theme(legend.text = element_text(size = 30), legend.title = element_text(size = 30))  +
  theme(legend.justification = c(1, 1), legend.position = c(1, 1), legend.background = element_rect(fill = NA)) + 
  theme(legend.key = element_rect(color = NA, fill = NA), legend.key.size = unit(0.5, 'inch')) +
  # modify the position of title
  theme(plot.title = element_text(hjust = 0.5, size = 50)) + 
  # modify the font size
  # modify the margin
  # theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
  theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), 'inch')) +
  theme(axis.text=element_text(size = 30), axis.title = element_text(size = 35), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1), axis.ticks.length = unit(0.12, 'inch')) 

dev.off()

mean(current_data_meta$cue[current_data_meta$cue_method == 1])
var(current_data_meta$cue[current_data_meta$cue_method == 1])


mean(current_data_meta$cue[current_data_meta$cue_method == 3])
var(current_data_meta$cue[current_data_meta$cue_method == 3])

#################################################################################
# PRODA results
#################################################################################
model_name = 'cesm2_clm5_mic_vr_v22'
nn_exp_name = 'exp_pc_cesm2_23'
time_domain = 'whole_time'

#-----------------------------------------------
# Load Projected SOC site by site
#-----------------------------------------------

para_names =  c('bio', 'cryo', 'q10', 'w_scaling', 
                'tau4cwd', 'tau4l1', 'tau4l2', 'tau4s1', 'tau4s2_death', 'tau4s2_enz', 'tau4s3', 'tau4s4', 
                'mm_const_assim', 'mm_const_decom', 
                'fcwdl2', 
                'pl1s1', 'pl2s1', 'pl3s4', 'l1_cue', 'l2l3_cue', 
                'mic_cue', 'pdeath2soc',
                'beta')


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

ss_para_result = readMat(paste(data_dir_output, 'mcmc_summary_', model_name, '/', model_name, '_para_mean.mat', sep = ''))
ss_para_result = ss_para_result$para.mean
colnames(ss_para_result) = para_names

env_info_ss = readMat(paste(data_dir_input, 'wosis_2019_snap_shot/wosis_2019_snapshot_hugelius_mishra_env_info.mat', sep = ''))
env_info_ss = env_info_ss$EnvInfo[ , 4:80]
colnames(env_info_ss) = grid_var_names

# site by site soc stock
da_summary = readMat(paste(data_dir_output, 'mcmc_summary_', model_name, '/', model_name, '_da_summary_mic.mat', sep = ''))
da_summary = da_summary$da.summary.mic
colnames(da_summary) = c('profile_num', 'mic_cue', 
                         'soc_0_30cm', 'soc_30_100cm', 'soc_100_200cm', 'mic_0_30cm', 
                         'mic_stock', 'enz_stock', 'doc_stock', 'poc_stock', 'soc_total', 
                         'bulk_A', 'bulk_I', 'bulk_K', 'bulk_V', 'bulk_Xi', 'bulk_E', 'bulk_NPP')

soc_content_0_30cm = da_summary[ , 'soc_0_30cm']/0.3/env_info_ss[ , 'Bulk_Density_0cm'] # unit gC/kg, unit of bulk density:# unit kg/m3 
mic_content_0_30cm = da_summary[ , 'mic_0_30cm']/0.3/env_info_ss[ , 'Bulk_Density_0cm'] # unit gC/kg, unit of bulk density:# unit kg/m3 

# load climate info
global_climate = env_info_ss[ , 'Koppen_Climate_2018']
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


# load ESA land cover
var_landcover = env_info_ss[ , 'ESA_Land_Cover']

# load soil texture
var_texture = env_info_ss[ , 'Texture_USDA_0cm']

#-----------------------------------------------
# all valid points and corresponding CUE-SOC relationship
#-----------------------------------------------
valid_loc = which(is.na(da_summary[ , 'profile_num']) == 0)

current_data_proda_total = data.frame(cbind(valid_loc, 
                                            var_climate[valid_loc],
                                            da_summary[valid_loc, 'bulk_A'], 
                                            da_summary[valid_loc, 'mic_cue'], 
                                            soc_content_0_30cm[valid_loc], 
                                            mic_content_0_30cm[valid_loc], 
                                            env_info_ss[valid_loc, 'Bulk_Density_0cm'],
                                            env_info_ss[valid_loc, 'CEC_0cm'],
                                            env_info_ss[valid_loc, 'Clay_Content_0cm'],
                                            env_info_ss[valid_loc, 'cesm2_npp'],
                                            env_info_ss[valid_loc, 'Annual Mean Temperature'] 
))

colnames(current_data_proda_total) = c('profile_num', 'climate', 'cue', 'mic_cue', 'soc', 'mic', 'bd', 'cec', 'clay', 'npp', 'mat')

current_data_proda_total = current_data_proda_total[is.na(apply(current_data_proda_total, 1, sum, na.rm = FALSE)) == 0, ]

#--------------------------- mix model microbial vs non-microbial biomass
mix_model = lmer((mic) ~ mic_cue + (1 | climate), data = current_data_proda_total)
mix_model_summary_mic = summary(mix_model)
mix_model_summary_mic

r_squared_mixed = summary(lm(current_data_proda_total$mic ~ predict(mix_model)))
r_squared_mixed 


mix_model = lmer((soc-mic) ~ mic_cue + (1 | climate), data = current_data_proda_total)
mix_model_summary_mic = summary(mix_model)
mix_model_summary_mic

r_squared_mixed = summary(lm((current_data_proda_total$soc-current_data_proda_total$mic) ~ predict(mix_model)))
r_squared_mixed 
#--------------------------- mix model bulk density
mix_model = lmer(log10(soc) ~ mic_cue + log10(bd) + (mic_cue | climate), data = current_data_proda_total)
mix_model_summary_mic = summary(mix_model)
mix_model_summary_mic

r_squared_mixed = summary(lm((current_data_proda_total$soc-current_data_proda_total$mic) ~ predict(mix_model)))
r_squared_mixed 

#--------------------------- mix model CEC
mix_model = lmer(log10(soc) ~ mic_cue + (cec) + (mic_cue | climate), data = current_data_proda_total)
mix_model_summary_mic = summary(mix_model)
mix_model_summary_mic

r_squared_mixed = summary(lm((current_data_proda_total$soc-current_data_proda_total$mic) ~ predict(mix_model)))
r_squared_mixed 


#--------------------------- mix model clay content
mix_model = lmer(log10(soc) ~ mic_cue + (clay) + (mic_cue | climate), data = current_data_proda_total)
mix_model_summary_mic = summary(mix_model)
mix_model_summary_mic

r_squared_mixed = summary(lm((current_data_proda_total$soc-current_data_proda_total$mic) ~ predict(mix_model)))
r_squared_mixed 

#--------------------------- mix model NPP
mix_model = lmer(log10(soc) ~ mic_cue + log10(npp+1) + (mic_cue | climate), data = current_data_proda_total)
mix_model_summary_mic = summary(mix_model)
mix_model_summary_mic

r_squared_mixed = summary(lm((current_data_proda_total$soc-current_data_proda_total$mic) ~ predict(mix_model)))
r_squared_mixed 

#--------------------------- final mixed model 1
mix_model = lmer(log10(soc) ~ mic_cue + (1 | climate), data = current_data_proda_total)
mix_model_summary_proda_total = summary(mix_model)
mix_model_summary_proda_total

r_squared_mixed = summary(lm(log10(current_data_proda_total$soc) ~ predict(mix_model)))
r_squared_mixed 

fix_effect_model = summary(mix_model_summary_proda_total)
fix_effect_pred = fix_effect_model$coefficients[1, 1] + 
  fix_effect_model$coefficients[2, 1]*current_data_proda_total$mic_cue
r_squared_fixed = summary(lm(log10(current_data_proda_total$soc) ~ fix_effect_pred))
r_squared_fixed
#--------------------------- final mixed model 2
mix_model = lmer(log10(soc) ~ mic_cue + (mic_cue | climate), data = current_data_proda_total)
mix_model_summary_proda_total = summary(mix_model)
mix_model_summary_proda_total
fit_function_proda_total = function(x) {10**(mix_model_summary_proda_total$coefficients[1, 1] + x*(mix_model_summary_proda_total$coefficients[2, 1]))}

r_squared_mixed = summary(lm(log10(current_data_proda_total$soc) ~ predict(mix_model)))
r_squared_mixed 

fix_effect_model = summary(mix_model_summary_proda_total)
fix_effect_pred = fix_effect_model$coefficients[1, 1] + 
  fix_effect_model$coefficients[2, 1]*current_data_proda_total$mic_cue
r_squared_fixed = summary(lm(log10(current_data_proda_total$soc) ~ fix_effect_pred))
r_squared_fixed

if (mix_model_summary_proda_total$coefficients[1, 5] < 0.001) {
  p_intercept = paste(' < 0.001', sep = '')
} else {
  p_intercept = paste(' = ', round(mix_model_summary_proda_total$coefficients[1, 5], 3), sep = '')
}

if (mix_model_summary_proda_total$coefficients[2, 5] < 0.001) {
  p_slope = paste(' < 0.001', sep = '')
} else {
  p_slope = paste(' = ', round(mix_model_summary_proda_total$coefficients[2, 5], 3), sep = '')
}

text_data = data.frame('x_axis' = c(0.0), 
                       'y_axis' = c(0.5), 
                       'equation' = paste('log10(SOC) = ', round(mix_model_summary_proda_total$coefficients[1, 1], 2), ' + ', round(mix_model_summary_proda_total$coefficients[2, 1], 2), '*CUE', '\nP(Intercept)', p_intercept, ', P(CUE)', p_slope, '\nExplained Variation = ', round(r_squared_mixed$r.squared*100, 2), '%', sep = ''))

#---------------------------- Plot Figures SOC content and cue (0 - 30cm)
p_proda_cue_soc_total =
  ggplot(data = current_data_proda_total) + 
  stat_bin_hex(aes(x = mic_cue, y = soc), bins = 100) +
  scale_fill_gradientn(name = 'Count', colors = viridis(7), trans = 'identity', limits = c(1, 50), oob = scales::squish) +
  scale_x_continuous(limits = c(0, NA), trans = 'identity') +
  scale_y_continuous(limits = c(0.1, 1000), trans = 'log10', labels = trans_format('log10', math_format(10^.x))) + 
  geom_function(fun = fit_function_proda_total, size = 2, color = 'black') + 
  geom_text(data = text_data, aes(x = x_axis, y = y_axis, label = equation), vjust = 1, hjust = 0, size = 9, show.legend = FALSE) + 
  theme_classic() + 
  # add title
  labs(title = '', x = 'Carbon use efficiency', y = expression(paste('Soil organic carbon (g C kg'^'-1', ')', sep = ''))) + 
  # change the legend properties
  guides(fill = guide_colorbar(direction = 'horizontal', barwidth = 15, barheight = 2.5, title.position = 'right', title.hjust = 0, title.vjust = 0.8, label.hjust = 0.5, frame.linewidth = 0), reverse = FALSE) +
  theme(legend.text = element_text(size = 25), legend.title = element_text(size = 25))  +
  theme(legend.justification = c(0, 1), legend.position = c(0, 1), legend.background = element_rect(fill = NA)) + 
  # modify the position of title
  theme(plot.title = element_text(hjust = 0.5, size = 50)) + 
  # modify the font size
  # modify the margin
  # theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
  theme(plot.margin = unit(c(0., 0.2, 0.2, 0.2), 'inch')) +
  theme(axis.text=element_text(size = 30), axis.title = element_text(size = 35), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1), axis.ticks.length = unit(0.12, 'inch')) 



#----------------------- Plot Figures MAT-CUE
corr_test = cor.test(current_data_proda_total$mat, current_data_proda_total$mic_cue)

if (corr_test$p.value < 0.001) { 
  text_data = data.frame('x_axis' = c(-15), 
                         'y_axis' = c(0.7), 
                         'equation' = paste('Pearson corelation = ', round(corr_test$estimate, 2), 
                                            '\nP-value < 0.001', 
                                            sep = ''))
} else {
  text_data = data.frame('x_axis' = c(-15), 
                         'y_axis' = c(0.7), 
                         'equation' = paste('Pearson corelation = ', round(corr_test$estimate, 2), 
                                            '\nP-value = ', round(corr_test$p.value, 3), 
                                            sep = ''))
}

p_proda_mat_cue_total =
  ggplot() +
    # stat_bin_hex(data = current_data_proda_total, aes(x = mat, y = mic_cue), bins = 100) +
    geom_point(data = current_data_proda_total, aes(x = mat, y = mic_cue), shape = 16, size = 0.5, color = 'snow4') +
    geom_smooth(data = current_data_proda_total, aes(x = mat, y = mic_cue), size = 2, color = 'black', fill = 'black', method = 'lm', alpha = 0.2) +
  scale_fill_gradientn(name = 'Count', colors = viridis(7), trans = 'identity', limits = c(1, 50), oob = scales::squish) +
  scale_x_continuous(trans = 'identity') +
  # change the background to black and white
  geom_text(data = text_data, aes(x = x_axis, y = y_axis, label = equation), vjust = 1, hjust = 0, size = 9, show.legend = FALSE) + 
  theme_classic() + 
  # add title
  labs(title = '', y = 'Carbon use efficiency', x = 'Mean annual temperature (°C)') + 
  # change the legend properties
  guides(fill = guide_colorbar(direction = 'horizontal', barwidth = 15, barheight = 2.5, title.position = 'right', title.hjust = 0, title.vjust = 0.8, label.hjust = 0.5, frame.linewidth = 0), reverse = FALSE) +
  theme(legend.text = element_text(size = 25), legend.title = element_text(size = 25))  +
  theme(legend.justification = c(0, 0), legend.position = c(0, 0), legend.background = element_rect(fill = NA)) + 
  # modify the position of title
  theme(plot.title = element_text(hjust = 0.5, size = 50)) + 
  # modify the font size
  # modify the margin
  # theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
  theme(plot.margin = unit(c(0., 0.2, 0.2, 0.2), 'inch')) +
  theme(axis.text=element_text(size = 30), axis.title = element_text(size = 35), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1), axis.ticks.length = unit(0.12, 'inch')) 


#-------------CUE at different scale relationship
# lm_model = lm(cue ~ mic_cue, data = current_data_proda_total)
# lm_model_summary_proda_total = summary(lm_model)
# lm_model_summary_proda_total
# 
# text_data = data.frame("x_axis" = c(0.01, 0.01), 
#                        "y_axis" = c(1, 0.92)*(0.5 - 0.1) + 0.1, 
#                        "equation" = c(paste("CUE['system']~' = '~", round(lm_model_summary_proda_total$coefficients[1, 1], 2), " + ", round(lm_model_summary_proda_total$coefficients[2, 1], 2), "~'*'~eta['DOC']", sep = ""), paste("Explained~~Variation~' = '~", round(lm_model_summary_proda_total$r.squared*100, 2), "~'%'", sep = "")))

diff_cue_corr = cor.test(current_data_proda_total$cue, current_data_proda_total$mic_cue)

text_data = data.frame('x_axis' = c(0.01), 
                       'y_axis' = c(0.6), 
                       'equation' = paste('Pearson corelation = ', round(diff_cue_corr$estimate, 2), 
                                          '\nP-value < 0.001', 
                                          sep = ''))

p_proda_total_diff_cue =
  ggplot(data = current_data_proda_total) + 
  stat_bin_hex(aes(x = mic_cue, y = cue), bins = 100) +
  scale_fill_gradientn(name = 'Count', colors = viridis(7), trans = 'identity', limits = c(1, 50), oob = scales::squish) +
  geom_smooth(aes(x = mic_cue, y = cue), size = 2, color = 'black', fill = 'black', method = 'lm', alpha = 0.2) +
  scale_x_continuous(trans = 'identity', n.breaks = 7, limits = c(0.0, 0.7)) +
  scale_y_continuous(trans = 'identity', n.breaks = 7, limits = c(NA, NA)) +
  # change the background to black and white
  geom_text(data = text_data, aes(x = x_axis, y = y_axis, label = equation), parse = FALSE, vjust = 1, hjust = 0, size = 9, show.legend = FALSE) + 
  # change the background to black and white
  theme_classic() +
  # add title
  labs(x = expression(paste(eta['DOC'], sep = '')), y = 'System level CUE') +
  # modify the position of title
  theme(plot.title = element_text(hjust = 0.5, size = 50)) + 
  theme(legend.justification = c(1, 0), legend.position = c(1, 0), legend.background = element_rect(fill = NA), legend.text.align = 0) +
  guides(fill = guide_colorbar(direction = 'horizontal', barwidth = 15, barheight = 2.5, title.position = 'right', title.hjust = 0, title.vjust = 0.8, label.hjust = 0.5, frame.linewidth = 0), reverse = FALSE) +
  theme(legend.text = element_text(size = 30), legend.title = element_text(size = 30))  +
  theme(legend.key = element_rect(color = NA, fill = NA), legend.key.size = unit(0.5, 'inch')) +
  theme(legend.direction = 'horizontal') + 
  # modify the margin
  # theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
  theme(plot.margin = unit(c(0., 0.2, 0.2, 0.2), 'inch')) +
  theme(axis.text=element_text(size = 30), axis.title = element_text(size = 35), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1), axis.ticks.length = unit(0.12, 'inch')) 


#################################################################################
# final figure
#################################################################################
# CUE-SOC
jpeg(paste('./Ensemble/revision4_main_fig2.jpeg', sep = ''), width = 18, height = 9, units = 'in', res = 300)

p_proda_cue_soc_total = p_proda_cue_soc_total + 
  labs(y = '  ')

plot_grid(p_meta_cue_soc, p_proda_cue_soc_total,
          labels = c('a', 'b'),
          label_size = 50,
          label_x = 0.02, label_y = 1.0,
          label_fontfamily = 'Arial',
          label_fontface = 'bold',
          nrow = 1)
dev.off()


pdf(paste('./Ensemble/revision4_main_fig2.pdf', sep = ''), width = 18, height = 9)

p_proda_cue_soc_total = p_proda_cue_soc_total + 
  labs(y = '  ')

plot_grid(p_meta_cue_soc, p_proda_cue_soc_total,
          labels = c('a', 'b'),
          label_size = 50,
          label_x = 0.02, label_y = 1.0,
          label_fontface = 'bold',
          nrow = 1)
dev.off()

# MAT-CUE
jpeg(paste('./Ensemble/revision4_mat_cue.jpeg', sep = ''), width = 18, height = 9, units = 'in', res = 300)

p_proda_mat_cue_total = p_proda_mat_cue_total + 
  labs(y = '  ')

plot_grid(p_meta_mat_cue, p_proda_mat_cue_total,
          labels = c('a', 'b'),
          label_size = 50,
          label_x = 0.02, label_y = 1.0,
          label_fontfamily = 'Arial',
          label_fontface = 'bold',
          nrow = 1)
dev.off()

# different scales of CUE

jpeg(paste('./Ensemble/revision4_diff_cue_mic.jpeg', sep = ''), width = 10, height = 10, units = 'in', res = 300)
p_proda_total_diff_cue
dev.off()

