## Packages
library(R.matlab)
library(maps)
library(ggplot2)
library(interp)
library(gridExtra)
library(grid)
library(ModelMetrics)
library(R.matlab)
library(jcolors)
library(Rfast)
library(scales)

library(cowplot)

library(ggcorrplot)
library(corrplot)

# library(colorplaner)

# dev.off()
##
rm(list = ls())

setwd('/Users/phoenix/Google_Drive/R')

jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
diff.colors <- colorRampPalette(c("#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B"))


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
# Data Pathway
#################################################################################
model_name = 'cesm2_clm5_mic_vr_v22'
nn_exp_name = 'exp_pc_cesm2_23' 
time_domain = 'whole_time'

data_dir_output = '/Users/phoenix/Google_Drive/Tsinghua_Luo/Projects/DATAHUB/ENSEMBLE/OUTPUT_DATA/'
data_dir_input = '/Users/phoenix/Google_Drive/Tsinghua_Luo/Projects/DATAHUB/ENSEMBLE/INPUT_DATA/'


#################################################################################
# Importance Cate to Cate
#################################################################################

var_category = c('Soil Structure', 'Soil Chemical', 'Climate', 'Vegetation', 'Geography')
para_category = c('Microbial CUE', 'Baseline Decomposition', 'Environmental Modifier', 'Vertical Transport', 'Input Allocation', 'Carbon Transfer Fraction')

ipara = 2
for (ipara in 1:length(para_category)) {
  mse_summary = read.csv(paste(data_dir_output, 'neural_networking/permutation_test_mse_summary_cate2cate_best_guess_var_', (ipara-1), '_', model_name, '_', time_domain, '_', nn_exp_name, '_cross_valid_0_7.csv', sep = ''), header = FALSE)
  
  relative_importance = mse_summary

    # for (ivar in 1:length(var_category)) {
    #   relative_importance[ , ivar] = mse_summary[ , ivar]/apply(mse_summary, 1, sum, na.rm = TRUE)
    # }
  
  current_data = data.frame(array(NA, dim = c(length(var_category), 5)))
  colnames(current_data) = c('Importance', 'Parameter', 'Category', 'Upper', 'Lower')
  
  for (ivar in 1:length(var_category)) {
    current_data[ivar, 'Importance'] = median(relative_importance[ , ivar], na.rm = TRUE)
    current_data[ivar, 'Lower'] = quantile(relative_importance[ , ivar], na.rm = TRUE, probs = 0.16) # mean(relative_importance[ , ivar], na.rm = TRUE) - sd(relative_importance[ , ivar], na.rm = TRUE) # quantile(relative_importance[ , ivar], na.rm = TRUE, probs = 0.025)
    current_data[ivar, 'Upper'] = quantile(relative_importance[ , ivar], na.rm = TRUE, probs = 0.84) # mean(relative_importance[ , ivar], na.rm = TRUE) + sd(relative_importance[ , ivar], na.rm = TRUE) #quantile(relative_importance[ , ivar], na.rm = TRUE, probs = 0.975)
    current_data[ivar, 'Parameter'] = para_category[ipara]
    current_data[ivar, 'Category'] = var_category[ivar]
    
  }
  
  
  sub_current_data = current_data[current_data$Parameter == para_category[ipara], ]
  
  p =
    ggplot(data = current_data, aes(x = reorder(Category, Importance), y = Importance)) + 
    geom_bar(width = 0.8, alpha = 1, linetype = 'solid', size = 0.7, color = 'black', fill = 'blue4', stat = 'identity', position = position_dodge(0.8)) +
    # geom_point(shape = 21, size = 7, color = 'white', stat = 'identity', position = position_dodge(0.8)) +
    geom_errorbar(aes(x = reorder(Category, Importance), ymin = Lower, ymax = Upper), width = 0.3, size = 0.7, stat = 'identity', position = position_dodge(0.8)) + 
    coord_flip() +
    theme_classic() +
    scale_x_discrete(name = '') +
    scale_y_continuous(limits=c(1, NA), oob = rescale_none) +
    theme(axis.text.x = element_text(hjust = 0.5, vjust = 0.5)) + 
    theme(axis.text.y = element_text(hjust = 0.5, vjust = 0.5)) +
    # add title
    labs(title = para_category[ipara], x = '', y = 'Permutation Importance') + 
    # modify the position of title
    theme(plot.title = element_text(hjust = 0.5, size = 25)) + 
    # modify the margin
    theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), 'inch')) +
    theme(axis.text=element_text(size = 20), axis.title = element_text(size = 20), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1), axis.ticks.length = unit(0.1, 'inch')) 
  
  
  eval(parse(text = paste('p', ipara, ' = p', sep = '')))
  
}


jpeg(paste('./Ensemble/revision4_relative_importance_env_var_to_process.jpeg', sep = ''), width = 18, height = 10, units = 'in', res = 300)
plot_grid(p1, p2, p3, p4, p5, p6,
          nrow = 2, ncol = 3)
dev.off()

