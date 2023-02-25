
library(ncdf4)
library(R.matlab)
library(ggplot2)
# library(qlcMatrix)
library(jcolors)
library(cowplot)

# dev.off()
rm(list = ls())

setwd('/Users/phoenix/Google_Drive/R/Ensemble')

model_name = 'cesm2_clm5_mic_vr_v22'

ncfname = '/Users/phoenix/Google_Drive/Tsinghua_Luo/Projects/DATAHUB/ENSEMBLE/INPUT_DATA/wosis_2019_snap_shot/soc_profile_wosis_2019_snapshot_hugelius_mishra.nc'
profile_info = nc_open(ncfname)
profile_info = ncvar_get(profile_info)

max_depth = readMat('/Users/phoenix/Google_Drive/Tsinghua_Luo/Projects/DATAHUB/ENSEMBLE/INPUT_DATA/wosis_2019_snap_shot/wosis_2019_snapshot_hugelius_mishra_env_info.mat')
max_depth = max_depth$EnvInfo[ , 73]

mcmc_results_r2 = readMat(paste('/Users/phoenix/Google_Drive/Tsinghua_Luo/Projects/DATAHUB/ENSEMBLE/OUTPUT_DATA/mcmc_summary_', model_name, '/', model_name, '_stat_r2.mat', sep = ''))
mcmc_results_r2 = mcmc_results_r2$stat.r2

mcmc_results_gr = readMat(paste('/Users/phoenix/Google_Drive/Tsinghua_Luo/Projects/DATAHUB/ENSEMBLE/OUTPUT_DATA/mcmc_summary_', model_name, '/', model_name, '_para_gr.mat', sep = ''))
mcmc_results_gr = mcmc_results_gr$para.gr

mcmc_results_std = readMat(paste('/Users/phoenix/Google_Drive/Tsinghua_Luo/Projects/DATAHUB/ENSEMBLE/OUTPUT_DATA/mcmc_summary_', model_name, '/', model_name, '_para_var.mat', sep = ''))
mcmc_results_std = mcmc_results_std$para.var


current_data = data.frame(cbind(profile_info[ , c(1, 4, 5, 6)], apply(mcmc_results_r2, 1, max, na.rm = TRUE), apply(mcmc_results_gr, 1, mean, na.rm = TRUE), max_depth, NA))
colnames(current_data) = c('id', 'lon', 'lat', 'layer_num', 'r2', 'gr', 'max_depth', 'source')


## Jet colorbar function
jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
# can be changed to state or world to have US and world map
world_coastline = rgdal::readOGR(dsn='/Users/phoenix/Google_Drive/Tsinghua_Luo/World_Vector_Shape/ne110m/ne_110m_land.shp',layer = 'ne_110m_land')
world_coastline <- fortify(world_coastline)
Map.Using = world_coastline


valid_site_wosis = which(current_data$id > 10000 & current_data$r2 > 0.0 & current_data$max_depth > 50 & current_data$gr < 1.05)
valid_site_mishra = which(current_data$id < 10000 & current_data$r2 > 0.0 & current_data$max_depth > 50 & current_data$gr < 1.05)
invalid_site = which(current_data$r2 <= 0.0 | current_data$max_depth <= 50 | current_data$gr >= 1.05)

current_data$source[valid_site_wosis] = 'A_WoSIS'
current_data$source[valid_site_mishra] = 'B_Mishra & Hugelius 2020'
current_data$source[invalid_site] = 'C_Invalid'


current_data = current_data[order(current_data$source, decreasing = TRUE), ]

p_wosis = 
  ggplot(data = current_data) +
  geom_point(aes(x = lon, y = lat, color = source), shape = 16, size = 0.2, alpha = 1) + 
  scale_color_manual(name = '', values = c('#005AB5', '#DC3220', 'grey'), labels = c('WoSIS 2019 Snapshot', 'Mishra & Hugelius 2020', 'Invalid')) +
  geom_polygon(data = Map.Using, aes(x = long, y = lat, group = group), fill = NA, color = 'black', size = 0.3) +
  ylim(c(-56, 80)) +
  # change the background to black and white
  theme_bw() +
  # change the legend properties
  # theme(legend.position = 'none') +
  theme(legend.justification = c(0, 0), legend.position = c(0, 0), legend.background = element_rect(fill = NA), legend.text.align = 0) +
  theme(legend.text = element_text(size = 35), legend.title = element_text(size = 35))  +
  guides(colour = guide_legend(override.aes = list(size = 5))) +
  theme(legend.text = element_text(size = 15), legend.title = element_text(size = 20)) +
  # add title
  labs(x = '', y = '') + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
  # modify the position of title
  theme(plot.title = element_text(hjust = 0.5, size = 40)) + 
  # modify the font size
  theme(axis.title = element_text(size = 20)) + 
  # modify the margin
  theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), 'inch')) +
  theme(axis.text=element_text(size = 30))


##################################################
## load emperical data
##################################################
cue_meta_info = readMat('/Users/phoenix/Google_Drive/Tsinghua_Luo/Projects/DATAHUB/CUE_Synthesis/cue_meta.mat')
cue_meta_info = cue_meta_info$meta.site.info
colnames(cue_meta_info) = c('site_id', 'lon', 'lat', 'mat', 'incubation_temp', 'map', 'depth', 'cue', 'mic', 'soc', 'pH', 'cue_method', 'source_id')

cue_meta_info[cue_meta_info[ , 'cue_method'] < 3, 'cue_method'] = 1

valid_loc = which(is.na(cue_meta_info[ , 'cue']) == 0 & is.na(cue_meta_info[ , 'soc']) == 0 & cue_meta_info[ , 'source_id'] < 1000)

####################################################################
# Geographical Location
####################################################################
world_coastline = rgdal::readOGR(dsn='/Users/phoenix/Google_Drive/Tsinghua_Luo/World_Vector_Shape/ne110m/ne_110m_land.shp',layer = 'ne_110m_land')
world_coastline <- fortify(world_coastline)
Map.Using = world_coastline

current_data_middle = cue_meta_info[valid_loc, ]
site_list = unique(current_data_middle[ , 1])
current_data = data.frame(array(NA, dim = c(length(site_list), 3)))
colnames(current_data) = c('lon', 'lat', 'num')

isite = 1
for (isite in 1:length(site_list)) {
  site_loc = which(current_data_middle[ , 1] == site_list[isite])
  current_data[isite, 'lon'] = current_data_middle[site_loc[1], 'lon']
  current_data[isite, 'lat'] = current_data_middle[site_loc[1], 'lat']
  current_data[isite, 'num'] = length(site_loc)
}

p_meta = 
  ggplot() +
  geom_point(data = current_data, aes(x = lon, y = lat, size = num), color = 'red3', alpha = 0.7, na.rm = TRUE) +
  scale_size_continuous(name = 'Record Number', breaks = c(2, 5, 8)) + 
  geom_polygon(data = Map.Using, aes(x = long, y = lat, group = group), fill = NA, color = 'black', size = 0.5) +
  # change the background to black and white
  theme_bw() +
  ylim(-56, 80) +
  # change the legend properties
  # theme(legend.position = 'none') +
  theme(legend.justification = c(0, 0), legend.position = c(0, 0), legend.background = element_rect(fill = NA), legend.text.align = 0) +
  # theme(legend.justification = c(0.5, 0), legend.position = c(0.5, 0), legend.background = element_rect(fill = NA), legend.direction = 'horizontal') +
  # change the size of colorbar
  guides(fill = guide_colorbar(barwidth = 2.5, barheight = 14), reverse = FALSE) + 
  theme(legend.text = element_text(size = 20), legend.title = element_text(size = 20)) +
  # add title
  labs(x = '', y = '') + 
  # modify the position of title
  theme(plot.title = element_text(hjust = 0.5, size = 50)) + 
  # modify the font size
  theme(axis.title = element_text(size = 20)) + 
  # modify the margin
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
  theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), 'inch')) +
  theme(axis.text=element_text(size = 35))



jpeg(paste('./profile_loc.jpeg', sep = ''), width = 10, height = 8, units = 'in', res = 300)
plot_grid(p_meta, p_wosis,
          nrow = 2, ncol = 1,
          labels = c('a', 'b'),
          label_size = 30,
          label_x = 0.0, label_y = 1.02,
          label_fontfamily = 'Arial',
          label_fontface = 'bold'
)

dev.off()
