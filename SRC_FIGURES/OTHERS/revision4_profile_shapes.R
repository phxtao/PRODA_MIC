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
library(ncdf4)

# dev.off()
##
rm(list = ls())

setwd('/Users/phoenix/Google_Drive/R')

#############################################################################
# Data Path
#############################################################################
model_name = 'cesm2_clm5_mic_vr_v22'
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

#############################################################################
# data loading
#############################################################################
profile_info_colname = c('profile_id', 'country_id', 'country_name', 'lat', 'lon', 'layer_num', 'date')
soc_data_colname = c('profile_id', 'date', 'upper_depth', 'lower_depth', 'node_depth', 'soc_weight', 'soc_stock', 'bulk_denstiy', 'is_pedo')

# soc content data
ncfname = paste(data_dir_input, 'wosis_2019_snap_shot/soc_data_integrate_wosis_2019_snapshot_hugelius_mishra.nc', sep = '')
wosis_soc_info = nc_open(ncfname)
wosis_soc_info = ncvar_get(wosis_soc_info)
colnames(wosis_soc_info) = soc_data_colname
# soc profile info
ncfname = paste(data_dir_input, 'wosis_2019_snap_shot/soc_profile_wosis_2019_snapshot_hugelius_mishra.nc', sep = '')
wosis_profile_info = nc_open(ncfname)
wosis_profile_info = ncvar_get(wosis_profile_info)
colnames(wosis_profile_info) = profile_info_colname

# sampled profiles
eligible_loc = readMat(paste(data_dir_input, 'data4nn/eligible_profile_loc_0_cesm2_clm5_cen_vr_v2_whole_time.mat', sep = ''))
eligible_loc = eligible_loc$eligible.loc.0

stat_r2 = readMat(paste(data_dir_output, 'mcmc_summary_', model_name, '/', model_name, '_stat_r2.mat', sep = ''))
stat_r2 = stat_r2$stat.r2
stat_r2 = apply(stat_r2, 1, max, na.rm = TRUE)

para_gr = readMat(paste(data_dir_output, 'mcmc_summary_', model_name, '/', model_name, '_para_gr.mat', sep = ''))
para_gr = para_gr$para.gr
para_gr = apply(para_gr, 1, mean, na.rm = TRUE)

valid_loc = intersect(which(stat_r2 > 0 & para_gr < 1.05), eligible_loc)
invalid_loc = setdiff(eligible_loc, valid_loc)

sample_num = 1000
set.seed(1)
sample_profile_loc = cbind(sample(valid_loc, sample_num, replace = FALSE), 
                           sample(invalid_loc, sample_num, replace = FALSE))

#############################################################################
# data processing
#############################################################################
current_data = c()

iset = 1
iprofile = 1
for (iset in 1:2) {
  
  for (iprofile in 1:length(sample_profile_loc[ , iset])) {
    print(paste('processing set ', iset, ' profile ', iprofile, sep = ''))
    
    profile_id = wosis_profile_info[sample_profile_loc[iprofile, iset], 'profile_id']
    profile_loc = which(wosis_soc_info[ , 'profile_id'] == profile_id)
    
    obs_depth = wosis_soc_info[profile_loc, 'node_depth']
    obs_soc = wosis_soc_info[profile_loc, 'soc_stock']
    
    valid_layer_loc = which(is.na(obs_depth) == 0 & is.na(obs_soc) == 0)
    obs_depth = obs_depth[valid_layer_loc]
    obs_soc = obs_soc[valid_layer_loc]
    
    
    depth_order = order(obs_depth)
    obs_depth = obs_depth[depth_order]
    obs_soc = obs_soc[depth_order]
    
    shape_feature = obs_soc/c(obs_soc[1], obs_soc[1:length(obs_soc)-1])
    
    if (length(which(shape_feature <= 1)) == length(obs_soc)) {
      shape_type = 1
    } else if (which(obs_soc == max(obs_soc)) != 1) {
      shape_type = 2
    } else {
      shape_type = 3
    }
    
    current_data_middle = cbind(profile_id, obs_depth, obs_soc, obs_soc/obs_soc[1], shape_type, iset)
    current_data = rbind(current_data, current_data_middle)
  }
  
}


#############################################################################
# plot figure
#############################################################################
color_scheme = c('#1E88E5', '#D81B60', '#FFC107')


current_data = data.frame(current_data)
colnames(current_data) = c('id', 'depth', 'soc', 'soc_norm', 'shape_type', 'set')
current_data$id[current_data$shape_type == 2] = current_data$id[current_data$shape_type == 2] + 10**8
current_data$id[current_data$shape_type == 3] = current_data$id[current_data$shape_type == 3] + 10**6

valid_loc = which(current_data$set == 1)
p1 =
  ggplot(data = current_data[valid_loc, ]) + 
  geom_line(aes(x = depth, y = soc_norm, group = as.factor(id), color = as.factor(shape_type)), alpha = 0.5, size = 1) +
  scale_y_continuous(trans = 'identity') +
  scale_x_continuous(trans = 'identity', limits = c(0, 300)) +
  scale_color_manual(name = '',  labels = c('Monotonic decrease', 'Highest in the middle', 'Irregular zigzag'), values = color_scheme) +
  # change the background to black and white
  theme_classic() +
  # add title
  labs(title = '', x = 'Soil depth (cm)', y = expression(paste('SOC content (normalised)', sep = ''))) + 
  # modify the position of title
  theme(plot.title = element_text(hjust = 0.5, size = 40)) + 
  theme(legend.justification = c(1, 1), legend.position = c(1, 1), legend.background = element_rect(fill = NA), legend.text.align = 0) +
  theme(legend.text = element_text(size = 30), legend.title = element_text(size = 30))  +
  theme(legend.key = element_rect(color = NA, fill = NA), legend.key.size = unit(0.5, 'inch')) +
  # theme(legend.direction = 'horizontal') + 
  # modify the margin
  # theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
  theme(plot.margin = unit(c(0., 0.2, 0.2, 0.2), 'inch')) +
  theme(axis.text=element_text(size = 30), axis.title = element_text(size = 35), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1), axis.ticks.length = unit(0.12, 'inch')) 

jpeg(paste('./Ensemble/revision4_profile_shape.jpeg', sep = ''), width = 9, height = 9, units = 'in', res = 300)
p1
dev.off()





