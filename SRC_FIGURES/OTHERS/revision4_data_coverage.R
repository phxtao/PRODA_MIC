library(plotbiomes) # https://github.com/valentinitnelav/plotbiomes
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

# dev.off()
##
rm(list = ls())

setwd('/Users/phoenix/Google_Drive/R')



#############################################################################
# Data Path
#############################################################################
data_dir_output = '/Users/phoenix/Google_Drive/Tsinghua_Luo/Projects/DATAHUB/ENSEMBLE/OUTPUT_DATA/'
data_dir_input = '/Users/phoenix/Google_Drive/Tsinghua_Luo/Projects/DATAHUB/ENSEMBLE/INPUT_DATA/'



#################################################################################
# meta-analysis data 
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

env_info = readMat('/Users/phoenix/Google_Drive/Tsinghua_Luo/Projects/DATAHUB/CUE_Synthesis/cue_meta_envinfo.mat')
env_info = env_info$EnvInfo[ , 4:80]
colnames(env_info) = grid_var_names
cue_meta_info = readMat('/Users/phoenix/Google_Drive/Tsinghua_Luo/Projects/DATAHUB/CUE_Synthesis/cue_meta.mat')
cue_meta_info = cue_meta_info$meta.site.info
colnames(cue_meta_info) = c('site_id', 'lon', 'lat', 'mat', 'incubation_temp', 'map', 'depth', 'cue', 'mic', 'soc', 'pH', 'cue_method', 'source_id')

# load climate info
global_climate = env_info[ , 'Koppen_Climate_2018']
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


# load soil order
global_soilorder = env_info[ , 'USDA_Suborder']
var_soilorder = array(NA, dim = c(length(global_soilorder), 1))

var_soilorder[] = 113 # others
var_soilorder[which(global_soilorder >= 5 & global_soilorder <= 7)] = 101 # Gelisols
var_soilorder[which(global_soilorder >= 10 & global_soilorder <= 13)] = 102 # Histosols
var_soilorder[which(global_soilorder >= 15 & global_soilorder <= 19)] = 103 # Spodosols
var_soilorder[which(global_soilorder >= 20 & global_soilorder <= 27)] = 104 # Andisols
var_soilorder[which(global_soilorder >= 30 & global_soilorder <= 34)] = 105 # Oxisols
var_soilorder[which(global_soilorder >= 40 & global_soilorder <= 45)] = 106 # Vertisols
var_soilorder[which(global_soilorder >= 50 & global_soilorder <= 56)] = 107 # Aridisols
var_soilorder[which(global_soilorder >= 60 & global_soilorder <= 64)] = 108 # Ultisols
var_soilorder[which(global_soilorder >= 69 & global_soilorder <= 77)] = 109 # Mollisols
var_soilorder[which(global_soilorder >= 80 & global_soilorder <= 84)] = 110 # Alfisols
var_soilorder[which(global_soilorder >= 85 & global_soilorder <= 86)] = 111 # Inceptisols
var_soilorder[which(global_soilorder >= 89 & global_soilorder <= 94)] = 111 # Inceptisols
var_soilorder[which(global_soilorder >= 95 & global_soilorder <= 99)] = 112 # Entisols

# load ESA land cover
global_landcover = env_info[ , 'ESA_Land_Cover']
var_landcover = array(NA, dim = c(length(global_landcover), 1))

var_landcover[which(global_landcover >= 1 & global_landcover <= 2)] = 101 # agriculture
var_landcover[which(global_landcover >= 3 & global_landcover <= 4)] = 102 # mosaic agriculture
var_landcover[which(global_landcover >= 5 & global_landcover <= 6)] = 103 # broadleaved forest 
var_landcover[which(global_landcover >= 7 & global_landcover <= 8)] = 104 # needleleaved forest 
var_landcover[which(global_landcover >= 9 & global_landcover <= 9)] = 105 # mixed forest 
var_landcover[which(global_landcover >= 10 & global_landcover <= 11)] = 106 # Mosaic tree and shrub
var_landcover[which(global_landcover >= 12 & global_landcover <= 12)] = 107 # shrub
var_landcover[which(global_landcover >= 13 & global_landcover <= 13)] = 108 # grassland
var_landcover[which(global_landcover >= 14 & global_landcover <= 15)] = 109 # Lichen & mosses, sparse vegatation
var_landcover[which(global_landcover >= 16 & global_landcover <= 18)] = 110 # wetland
var_landcover[which(global_landcover >=19)] = 111 # urban and other

valid_loc = which(cue_meta_info[ , 'source_id'] < 1000)
current_data_meta = data.frame(cbind(cue_meta_info[valid_loc, c('mat', 'map')],
                                     var_climate[valid_loc],
                                     env_info[valid_loc, c('Texture_USDA_0cm')],
                                     var_soilorder[valid_loc],
                                     var_landcover[valid_loc]
))

colnames(current_data_meta) = c('mat', 'map', 'climate', 'texture', 'order', 'land_cover')


#################################################################################
# PRODA results
#################################################################################
model_name = 'cesm2_clm5_mic_vr_v22'
nn_exp_name = 'exp_pc_cesm2_23'
time_domain = 'whole_time'

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
valid_loc = which(is.na(da_summary[ , 'profile_num']) == 0)

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

# load soil order
global_soilorder = env_info_ss[ , 'USDA_Suborder']
var_soilorder = array(NA, dim = c(length(global_soilorder), 1))

var_soilorder[] = 113 # others
var_soilorder[which(global_soilorder >= 5 & global_soilorder <= 7)] = 101 # Gelisols
var_soilorder[which(global_soilorder >= 10 & global_soilorder <= 13)] = 102 # Histosols
var_soilorder[which(global_soilorder >= 15 & global_soilorder <= 19)] = 103 # Spodosols
var_soilorder[which(global_soilorder >= 20 & global_soilorder <= 27)] = 104 # Andisols
var_soilorder[which(global_soilorder >= 30 & global_soilorder <= 34)] = 105 # Oxisols
var_soilorder[which(global_soilorder >= 40 & global_soilorder <= 45)] = 106 # Vertisols
var_soilorder[which(global_soilorder >= 50 & global_soilorder <= 56)] = 107 # Aridisols
var_soilorder[which(global_soilorder >= 60 & global_soilorder <= 64)] = 108 # Ultisols
var_soilorder[which(global_soilorder >= 69 & global_soilorder <= 77)] = 109 # Mollisols
var_soilorder[which(global_soilorder >= 80 & global_soilorder <= 84)] = 110 # Alfisols
var_soilorder[which(global_soilorder >= 85 & global_soilorder <= 86)] = 111 # Inceptisols
var_soilorder[which(global_soilorder >= 89 & global_soilorder <= 94)] = 111 # Inceptisols
var_soilorder[which(global_soilorder >= 95 & global_soilorder <= 99)] = 112 # Entisols

# load ESA land cover
global_landcover = env_info_ss[ , 'ESA_Land_Cover']
var_landcover = array(NA, dim = c(length(global_landcover), 1))

var_landcover[which(global_landcover >= 1 & global_landcover <= 2)] = 101 # agriculture
var_landcover[which(global_landcover >= 3 & global_landcover <= 4)] = 102 # mosaic agriculture
var_landcover[which(global_landcover >= 5 & global_landcover <= 6)] = 103 # broadleaved forest 
var_landcover[which(global_landcover >= 7 & global_landcover <= 8)] = 104 # needleleaved forest 
var_landcover[which(global_landcover >= 9 & global_landcover <= 9)] = 105 # mixed forest 
var_landcover[which(global_landcover >= 10 & global_landcover <= 11)] = 106 # Mosaic tree and shrub
var_landcover[which(global_landcover >= 12 & global_landcover <= 12)] = 107 # shrub
var_landcover[which(global_landcover >= 13 & global_landcover <= 13)] = 108 # grassland
var_landcover[which(global_landcover >= 14 & global_landcover <= 15)] = 109 # Lichen & mosses, sparse vegatation
var_landcover[which(global_landcover >= 16 & global_landcover <= 18)] = 110 # wetland
var_landcover[which(global_landcover >=19)] = 111 # urban and other

current_data_proda = data.frame(env_info_ss[valid_loc, c('Annual Mean Temperature', 'Annual Precipitation')],
                                var_climate[valid_loc],
                                env_info_ss[valid_loc, c('Texture_USDA_0cm')],
                                var_soilorder[valid_loc],
                                var_landcover[valid_loc])
colnames(current_data_proda) = c('mat', 'map', 'climate', 'texture', 'order', 'land_cover')


#################################################################################
# Raw data
#################################################################################
eligible_loc = readMat(paste(data_dir_input, 'data4nn/eligible_profile_loc_0_cesm2_clm5_cen_vr_v2_whole_time.mat', sep = ''))
eligible_loc = eligible_loc$eligible.loc.0

current_data_raw = data.frame(env_info_ss[eligible_loc, c('Annual Mean Temperature', 'Annual Precipitation')],
                              var_climate[eligible_loc],
                              env_info_ss[eligible_loc, c('Texture_USDA_0cm')],
                              var_soilorder[eligible_loc],
                              var_landcover[eligible_loc])
colnames(current_data_raw) = c('mat', 'map', 'climate', 'texture', 'order', 'land_cover')

#################################################################################
# plot figure land cover
#################################################################################
land_cover_var_list = c('Agriculture', 'Mosaic agriculture', 'Broadleaved forest', 'Needleleaved forest', 
                        'Mixed forest', 'Mosaic tree & shrub', 'Shrub', 'Grassland', 'Sparse vegatation',
                        'Wetland', 'Urban & other')
climate_var_list = c('Af', 'Am', 'Aw', 'BW', 'BS', 'Cs', 'Cw', 'Cf', 'Ds', 'Dw', 'Df', 'E')
texture_var_list = c('Cl', 'SiCl', 'SaCl', 'ClLo', 'SiClLo', 'SaClLo', 'Lo', 'SiLo', 'SaLo', 'Si', 'LoSa', 'Sa')
soilorder_var_list = c('Gelisols', 'Histosols', 'Spodosols', 'Andisols', 'Oxisols', 'Vertisols', 'Aridisols', 'Ultisols', 'Mollisols', 'Alfisols', 'Inceptisols', 'Entisols', 'Others')


current_data_summary_climate = array(NA, dim = c(12, 4))
current_data_summary_land_cover = array(NA, dim = c(11, 4))
current_data_summary_texture = array(NA, dim = c(12, 4))
current_data_summary_soilorder = array(NA, dim = c(13, 4))

ivar = 1
for (ivar in 1:13) {
  # soil order
  soilorder_id = 100+ivar
  current_data_summary_soilorder[ivar, ] = c(soilorder_id,
                                             length(which(current_data_raw$order == soilorder_id))/length(current_data_raw$order),
                                             length(which(current_data_proda$order == soilorder_id))/length(current_data_proda$order),
                                           length(which(current_data_meta$order == soilorder_id))/length(current_data_meta$order)
  )
  
  if (ivar < 13) {
    # climate
    climate_id = 100+ivar
    current_data_summary_climate[ivar, ] = c(climate_id,
                                             length(which(current_data_raw$climate == climate_id))/length(current_data_raw$climate),
                                             length(which(current_data_proda$climate == climate_id))/length(current_data_proda$climate),
                                             length(which(current_data_meta$climate == climate_id))/length(current_data_meta$climate)
    )
    
    # soil texture
    texture_id = ivar
    current_data_summary_texture[ivar, ] = c(texture_id,
                                             length(which(current_data_raw$texture == texture_id))/length(current_data_raw$texture),
                                             length(which(current_data_proda$texture == texture_id))/length(current_data_proda$texture),
                                             length(which(current_data_meta$texture == texture_id))/length(current_data_meta$texture)
    )
  }
  
  # land cover
  if (ivar < 12) {
    land_cover_id = 100+ivar
    current_data_summary_land_cover[ivar, ] = c(land_cover_id,
                                                length(which(current_data_raw$land_cover == land_cover_id))/length(current_data_raw$land_cover),
                                                length(which(current_data_proda$land_cover == land_cover_id))/length(current_data_proda$land_cover),
                                                length(which(current_data_meta$land_cover == land_cover_id))/length(current_data_meta$land_cover)
    )
    
  } 
}

# current data climate
current_data_plot_climate = rbind(cbind(current_data_summary_climate[ , c(1, 2)], 
                                        1, # raw
                                        1), # climate 
                                  cbind(current_data_summary_climate[ , c(1, 3)], 
                                        2, # proda
                                        1), # climate 
                                  cbind(current_data_summary_climate[ , c(1, 4)], 
                                        3, # meta
                                        1) # climate
)

current_data_plot_climate = data.frame(current_data_plot_climate)
colnames(current_data_plot_climate) = c('subclass', 'percent', 'source', 'class')

# current data land cover
current_data_plot_land_cover = rbind(cbind(current_data_summary_land_cover[ , c(1, 2)], 
                                           1, # raw
                                           1), # climate 
                                     cbind(current_data_summary_land_cover[ , c(1, 3)], 
                                           2, # proda
                                           1), # climate 
                                     cbind(current_data_summary_land_cover[ , c(1, 4)], 
                                           3, # meta
                                           1) # climate
)

current_data_plot_land_cover = data.frame(current_data_plot_land_cover)
colnames(current_data_plot_land_cover) = c('subclass', 'percent', 'source', 'class')

# current data texture
current_data_plot_texture = rbind(cbind(current_data_summary_texture[ , c(1, 2)], 
                                        1, # raw
                                        1), # climate 
                                  cbind(current_data_summary_texture[ , c(1, 3)], 
                                        2, # proda
                                        1), # climate 
                                  cbind(current_data_summary_texture[ , c(1, 4)], 
                                        3, # meta
                                        1) # climate
)

current_data_plot_texture = data.frame(current_data_plot_texture)
colnames(current_data_plot_texture) = c('subclass', 'percent', 'source', 'class')

# current data soil order
current_data_plot_soilorder = rbind(cbind(current_data_summary_soilorder[ , c(1, 2)], 
                                          1, # raw
                                          1), # climate 
                                    cbind(current_data_summary_soilorder[ , c(1, 3)], 
                                          2, # proda
                                          1), # climate 
                                    cbind(current_data_summary_soilorder[ , c(1, 4)], 
                                          3, # meta
                                          1) # climate
)


current_data_plot_soilorder = data.frame(current_data_plot_soilorder)
colnames(current_data_plot_soilorder) = c('subclass', 'percent', 'source', 'class')



color_scheme = c('#88CCEE', '#CC6677', '#DDCC77', '#117733', '#332288', '#AA4499', '#44AA99', '#999933', '#882255', '#661100', '#6699CC', '#888888', 'black')

##############################climate 
current_data_plot_climate = current_data_plot_climate[current_data_plot_climate$source != 1, ]

p_climate =
  ggplot() + 
  geom_bar(data = current_data_plot_climate, aes(x = as.factor(source), y = percent*100, fill = as.factor(subclass)), color = 'white', stat = 'identity', position = 'stack', width = 0.7, size = 1) +
  scale_x_discrete(labels = c('PRODA', 'Meta')) +
  scale_fill_manual(name = '', labels = climate_var_list, values = color_scheme) + 
  
  theme_classic() + 
  # add title
  labs(title = 'Climate', x = '', y = 'Percentage (%)') + 
  # change the legend properties
  theme(legend.text = element_text(size = 25), legend.title = element_text(size = 25))  +
  theme(legend.justification = 'right', legend.position = 'right', legend.background = element_rect(fill = NA)) + 
  theme(legend.key = element_rect(color = NA, fill = NA), legend.key.size = unit(0.5, 'inch')) +
  # modify the position of title
  theme(plot.title = element_text(hjust = 0.5, size = 50)) + 
  # modify the font size
  # modify the margin
  # theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
  theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), 'inch')) +
  theme(axis.text=element_text(size = 30), axis.title = element_text(size = 35), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1), axis.ticks.length = unit(0.12, 'inch')) 


##############################texture 
current_data_plot_texture = current_data_plot_texture[current_data_plot_texture$source != 1, ]

p_texture = 
  ggplot() + 
  geom_bar(data = current_data_plot_texture, aes(x = as.factor(source), y = percent*100, fill = as.factor(subclass)), color = 'white', stat = 'identity', position = 'stack', width = 0.7, size = 1) +
  scale_x_discrete(labels = c('PRODA', 'Meta')) +
  scale_fill_manual(name = '', labels = texture_var_list, values = color_scheme) + 
  theme_classic() + 
  # add title
  labs(title = 'Soil texture', x = '', y = 'Percentage (%)') + 
  # change the legend properties
  theme(legend.text = element_text(size = 25), legend.title = element_text(size = 25))  +
  theme(legend.justification = 'right', legend.position = 'right', legend.background = element_rect(fill = NA)) + 
  theme(legend.key = element_rect(color = NA, fill = NA), legend.key.size = unit(0.5, 'inch')) +
  # modify the position of title
  theme(plot.title = element_text(hjust = 0.5, size = 50)) + 
  # modify the font size
  # modify the margin
  # theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
  theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), 'inch')) +
  theme(axis.text=element_text(size = 30), axis.title = element_text(size = 35), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1), axis.ticks.length = unit(0.12, 'inch')) 

##############################land cover 
current_data_plot_land_cover = current_data_plot_land_cover[current_data_plot_land_cover$source != 1, ]

p_land_cover = 
  ggplot() + 
  geom_bar(data = current_data_plot_land_cover, aes(x = as.factor(source), y = percent*100, fill = as.factor(subclass)), color = 'white', stat = 'identity', position = 'stack', width = 0.7, size = 1) +
  scale_x_discrete(labels = c('PRODA', 'Meta')) +
  scale_fill_manual(name = '', labels = land_cover_var_list, values = color_scheme) + 
  theme_classic() + 
  # add title
  labs(title = 'Land cover', x = '', y = 'Percentage (%)') + 
  # change the legend properties
  theme(legend.text = element_text(size = 25), legend.title = element_text(size = 25))  +
  theme(legend.justification = 'right', legend.position = 'right', legend.background = element_rect(fill = NA)) + 
  theme(legend.key = element_rect(color = NA, fill = NA), legend.key.size = unit(0.5, 'inch')) +
  # modify the position of title
  theme(plot.title = element_text(hjust = 0.5, size = 50)) + 
  # modify the font size
  # modify the margin
  # theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
  theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), 'inch')) +
  theme(axis.text=element_text(size = 30), axis.title = element_text(size = 35), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1), axis.ticks.length = unit(0.12, 'inch')) 

##############################soil order
current_data_plot_soilorder = current_data_plot_soilorder[current_data_plot_soilorder$source != 1, ]

p_soil_order =
  ggplot() + 
  geom_bar(data = current_data_plot_soilorder, aes(x = as.factor(source), y = percent*100, fill = as.factor(subclass)), color = 'white', stat = 'identity', position = 'stack', width = 0.7, size = 1) +
  scale_x_discrete(labels = c('PRODA', 'Meta')) +
  scale_fill_manual(name = '', labels = soilorder_var_list, values = color_scheme) + 
  theme_classic() + 
  # add title
  labs(title = 'Soil order', x = '', y = 'Percentage (%)') + 
  # change the legend properties
  theme(legend.text = element_text(size = 25), legend.title = element_text(size = 25))  +
  theme(legend.justification = 'right', legend.position = 'right', legend.background = element_rect(fill = NA)) + 
  theme(legend.key = element_rect(color = NA, fill = NA), legend.key.size = unit(0.5, 'inch')) +
  # modify the position of title
  theme(plot.title = element_text(hjust = 0.5, size = 50)) + 
  # modify the font size
  # modify the margin
  # theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
  theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), 'inch')) +
  theme(axis.text=element_text(size = 30), axis.title = element_text(size = 35), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1), axis.ticks.length = unit(0.12, 'inch')) 


jpeg(paste('./Ensemble/revision4_data_coverage_multi-var.jpeg', sep = ''), width = 35, height = 10, units = 'in', res = 300)

plot_grid(p_climate, p_texture, p_soil_order, p_land_cover,
          labels = c('a', 'b', 'c', 'd'),
          rel_widths = c(3, 3.1, 3.3, 3.8),
          label_size = 60,
          label_x = 0.02, label_y = 1.03,
          label_fontfamily = 'Arial',
          label_fontface = 'bold',
          nrow = 1)
dev.off()

