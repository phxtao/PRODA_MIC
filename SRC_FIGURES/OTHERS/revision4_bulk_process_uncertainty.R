## Packages
library(R.matlab)
library(ggplot2)
library(cowplot)
library(jcolors)
library(gridExtra)
library(viridis)

library(GGally)

library(rgdal)
library(raster)
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
model_name = 'cesm2_clm5_mic_vr_v22'
nn_exp_name = 'exp_pc_cesm2_23'
time_domain = 'whole_time'

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
# Env. Var Names
#################################################################################

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

valid_grid_loc = read.csv(paste(data_dir_output, 'neural_networking/valid_grid_loc_', model_name, '_', time_domain, '_', nn_exp_name, '_cross_valid_0_', as.character(1), '.csv', sep = ''), header = FALSE)
valid_grid_loc = valid_grid_loc$V1
grid_env_info = readMat(paste(data_dir_input, 'data4nn/world_grid_envinfo_present.mat', sep = ''))
grid_env_info = grid_env_info$EnvInfo
colnames(grid_env_info) = grid_var_names

cesm2_npp_std = grid_env_info[valid_grid_loc, 'cesm2_npp_std']

#################################################################################
# Load Projected SOC PRODA
#################################################################################
bootstrap_num = 200

bulk_process_list = c('A', 'I', 'K', 'V', 'Xi', 'NPP')

global_lat_lon = readMat(paste(data_dir_output, 'world_simulation_analyses/marginal_sensitivity_proda_NPP_', model_name, '_', nn_exp_name, '_cross_valid_0_', as.character(1), '.mat', sep = ''))
global_lat_lon = global_lat_lon$var.data.middle[[1]][[1]][ , 1:2]

bulk_process_std = readMat(paste(data_dir_output, 'world_simulation_analyses/uncertainty_process_bootstrap_', model_name, '.mat', sep = ''))
bulk_process_std = bulk_process_std$uncertainty.process
colnames(bulk_process_std) = c('soc', 'A', 'K', 'Xi', 'V', 'I', 'E', 'NPP')
bulk_process_std = bulk_process_std[ , c('soc', 'A', 'I', 'K', 'V', 'Xi', 'E', 'NPP')]

bulk_process_std[ , 'NPP'] = cesm2_npp_std
#################################################################################
# Plot Figures
#################################################################################

bulk_process_list = c('A', 'I', 'K', 'V', 'Xi', 'E', 'NPP')

process_scale_option = c('identity', 'identity', 'identity', 'identity', 'identity', 'identity', 'identity')
process_axis_label = c('Microbial CUE', 
                       'Input Allocation',
                       expression(paste('Baseline Decomposition (yr'^'-1', ')', sep = '')), 
                       expression(paste('Vertical Transport Rate (yr'^'-1', ')', sep = '')),
                       'Environmental Modifer',
                       expression(paste('Carbon Transfer Fraction', sep = '')), 
                       expression(paste('Input Carbon (g C yr'^'-1', ')', sep = '')))
process_name =  c('Microbial CUE', 
                  'Input Allocation', 
                  'Baseline Decomposition', 
                  'Vertical Transport Rate',
                  'Environmental Modifer', 
                  'Carbon Transfer Fraction',
                  'Input Carbon')
process_unit = c('unitless',
                 'unitless', 
                 expression(paste('yr'^'-1', sep = '')),
                 expression(paste('yr'^'-1', sep = '')),
                 'unitless', 
                 'unitless',
                 expression(paste('g C yr'^'-1', sep = '')))


legend_limit_lower_uncertain = c(0.,   0.035, 0.0,   0,     0.0,  0,   0) # apply(bulk_process_std, 2, quantile, prob = 0.005, na.rm = TRUE)
legend_limit_upper_uncertain = c(0.04, 0.075, 0.1, 0.04,  0.04, 0.02, 100) # apply(bulk_process_std, 2, quantile, prob = 0.995, na.rm = TRUE)

##################################soc stock and Residence Time map
world_coastline = rgdal::readOGR(dsn='/Users/phoenix/Google_Drive/Tsinghua_Luo/World_Vector_Shape/ne110m/ne_110m_land.shp',layer = 'ne_110m_land')
world_coastline <- spTransform(world_coastline, CRS('+proj=robin'))
world_coastline <- fortify(world_coastline)
Map.Using = world_coastline

ocean_left = cbind(rep(-180, 100), seq(from = 80, to = -56, by = -(80 + 56)/(100 -1)))
ocean_right = cbind(rep(180, 100), seq(from = -56, to = 80, by = (80 + 56)/(100 -1)))
ocean_top = cbind(seq(from = 180, to = -180, by = -(360)/(100 -1)), rep(80, 100))
ocean_bottom = cbind(seq(from = -180, to = 180, by = (360)/(100 -1)), rep(-56, 100))

world_ocean = rbind(ocean_left, ocean_bottom, ocean_right, ocean_top)
world_ocean = as.matrix(world_ocean)

world_ocean <- project(xy = world_ocean, proj = '+proj=robin') 

world_ocean = data.frame(world_ocean)
colnames(world_ocean) = c('long', 'lat')

color_scheme = c('#D55E00', '#E69F00', '#009E73', '#0072B2')

# soc map
CurrentData = data.frame(cbind(global_lat_lon, bulk_process_std[ , 'soc']))
colnames(CurrentData) = c('Lon', 'Lat', 'std_soc')

lon_lat_transfer = project(xy = as.matrix(CurrentData[ , c('Lon', 'Lat')]), proj = '+proj=robin') 
CurrentData[ , c('Lon', 'Lat')] = lon_lat_transfer

lat_limits = rbind(c(0, -56), c(0, 80))
lat_limits_robin = project(xy = as.matrix(lat_limits), proj = '+proj=robin') 

p_std_soc =
  ggplot() +
  geom_tile(data = CurrentData, aes(x = Lon, y = Lat, fill = std_soc/1000), height = 60000, width = 60000, na.rm = TRUE) +
  scale_fill_gradientn(name = expression(paste('kg C m'^'-2', sep = '')), colours = rev(viridis(15)), na.value="transparent", limits = c(1, 30), trans = 'log10', oob = scales::squish) +
  geom_polygon(data = Map.Using, aes(x = long, y = lat, group = group), fill = NA, color = 'black', size = 0.5) +
  geom_polygon(data = world_ocean, aes(x = long, y = lat), fill = NA, color = 'black', size = 1) +
  # change the background to black and white
  coord_equal() +
  # theme_map() +
  ylim(lat_limits_robin[ , 2]) +
  # change the legend properties
  # theme(legend.position = 'none') +
  theme(legend.justification = c(0.1, 0.36), legend.position = c(0.1, 0.36), legend.background = element_rect(fill = NA), legend.text.align = 0) +
  # theme(legend.justification = c(0.5, 0), legend.position = c(0.5, 0), legend.background = element_rect(fill = NA), legend.direction = 'horizontal') +
  # change the size of colorbar
  guides(fill = guide_colorbar(direction = 'vertical', barwidth = 2.5, barheight = 14, title.position = 'top', title.hjust = 0, label.hjust = 0.5, frame.linewidth = 0), reverse = FALSE) +
  theme(legend.text = element_text(size = 20, ), legend.title = element_text(size = 30)) +
  # add title
  labs(title = 'Uncertainty SOC Stock', x = '', y = '') + 
  # modify the position of title
  theme(plot.title = element_text(hjust = 0.5, vjust = -1, size = 50)) + 
  # modify the font size
  theme(axis.title = element_text(size = 20)) + 
  theme(panel.background = element_rect(fill = NA, colour = NA)) +
  # modify the margin
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
  theme(plot.margin = unit(c(0, 0, 0, 0), 'inch')) +
  theme(axis.text=element_text(size = 35))

################################## process map

ipara = 1
## Projected para 
for (ipara in 1:length(process_name)){
  
  CurrentData = data.frame(cbind(global_lat_lon, bulk_process_std[ , (ipara+1)]))
  colnames(CurrentData) = c('Lon', 'Lat', 'Project')
  
  lon_lat_transfer = project(xy = as.matrix(CurrentData[ , c('Lon', 'Lat')]), proj = '+proj=robin') 
  CurrentData[ , c('Lon', 'Lat')] = lon_lat_transfer
  
  
  p =
    ggplot() +
    geom_tile(data = CurrentData, aes(x = Lon, y = Lat, fill = Project), height = 60000, width = 60000, na.rm = TRUE) +
    scale_fill_gradientn(name = process_unit[ipara], colours = rev(viridis(15)), na.value="transparent", limits = c(legend_limit_lower_uncertain[ipara], legend_limit_upper_uncertain[ipara]), trans = process_scale_option[ipara], oob = scales::squish) +
    geom_polygon(data = Map.Using, aes(x = long, y = lat, group = group), fill = NA, color = 'black', size = 0.5) +
    geom_polygon(data = world_ocean, aes(x = long, y = lat), fill = NA, color = 'black', size = 1) +
    # change the background to black and white
    coord_equal() +
    # theme_map() +
    ylim(lat_limits_robin[ , 2]) +
    # change the legend properties
    # theme(legend.position = 'none') +
    theme(legend.justification = c(0.1, 0.36), legend.position = c(0.1, 0.36), legend.background = element_rect(fill = NA), legend.text.align = 0) +
    # theme(legend.justification = c(0.5, 0), legend.position = c(0.5, 0), legend.background = element_rect(fill = NA), legend.direction = 'horizontal') +
    # change the size of colorbar
    guides(fill = guide_colorbar(direction = 'vertical', barwidth = 2.5, barheight = 14, title.position = 'top', title.hjust = 0, label.hjust = 0.5, frame.linewidth = 0), reverse = FALSE) +
    theme(legend.text = element_text(size = 20, ), legend.title = element_text(size = 30)) +
    # add title
    labs(title = paste('Uncertainty ', process_name[ipara], sep = ''), x = '', y = '') + 
    # modify the position of title
    theme(plot.title = element_text(hjust = 0.5, vjust = -1, size = 50)) + 
    # modify the font size
    theme(axis.title = element_text(size = 20)) + 
    theme(panel.background = element_rect(fill = NA, colour = NA)) +
    # modify the margin
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
    theme(plot.margin = unit(c(0, 0, 0, 0), 'inch')) +
    theme(axis.text=element_text(size = 35))
  
  eval(parse(text = paste('p_std', ipara, ' = p', sep = '')))
  
}


jpeg(paste('./Ensemble/revision4_process_std.jpeg', sep = ''), width = 27, height = 25, units = 'in', res = 300)
plot_grid(p_std_soc, p_std1, 
          p_std6, p_std7,
          p_std2, p_std3,
          p_std5, p_std4, 
          labels = c('a', 'b', 
                     'c', 'd',
                     'e', 'f',
                     'g', 'h'),
          label_size = 70,
          label_x = 0.02, label_y = 1.03,
          label_fontfamily = 'Arial',
          label_fontface = 'bold',
          nrow = 4)
dev.off()
