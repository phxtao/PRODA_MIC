## Packages
library(R.matlab)
library(maps)
library(ggplot2)
library(R.matlab)
library(cowplot)
library(jcolors)
library(viridis)
library(cowplot)
library(FactoMineR)
library(factoextra)
library(ggrepel)
library(corrplot)
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
# Sensitivity component control
#################################################################################
# color_scheme = c('#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2', '#D55E00', '#CC79A7')
color_scheme = c('#e41a1c', '#984ea3', '#ff7f00', '#4daf4a', '#377eb8', '#a65628', 'grey3')
# color_scheme = c('#88CCEE', '#CC6677', '#DDCC77', '#117733', '#332288', '#AA4499', '#44AA99', '#999933', '#882255', '#661100', '#6699CC', '#888888')
# color_scheme = c('#7F3C8D', '#11A579', '#3969AC', '#F2B701', '#E73F74', '#80BA5A', '#E68310', '#008695', '#CF1C90', '#f97b72', '#4b4b8f', '#A5AA99')


bootstrap_num = 200
stat_matric = c('coeff_efficiency', 'correlation', 'coeff_determination', 'mse')

component_name = c('H_NN', 'A_A', 'F_I', 'D_K', 'E_V', 'B_Xi', 'C_NPP', 'A_E', 'Ga_Homo')

depth_name = c('C_30cm', 'B_100cm', 'A_200cm', 'D_full_depth')

control_summary = readMat(paste(data_dir_output, 'world_simulation_analyses/compontent_control_summary_bootstrap_', model_name, '.mat', sep = ''))
control_summary = control_summary$control.summary
control_summary = data.frame(control_summary)

colnames(control_summary) = c('var_soc', 'upper_soc', 'lower_soc', 'var_stat', 'upper_stat', 'lower_stat', 'depth', 'component')

for (icomponent in 1:length(component_name)) {
  control_summary[control_summary[ , 'component'] == icomponent, 'component'] = component_name[icomponent]
}

for (idepth in 1:length(depth_name)) {
  control_summary[control_summary[ , 'depth'] == idepth, 'depth'] = depth_name[idepth]
}

### plot figure component only
current_data_range = control_summary[which(control_summary$depth == 'B_100cm' & control_summary$component == 'H_NN'), ]
current_data = control_summary[which(control_summary$depth == 'B_100cm' & control_summary$component != 'Ga_Homo' & control_summary$component != 'H_NN'), ]

abs(current_data_range$var_soc - current_data$var_soc)
mean(abs(current_data_range$var_soc - current_data$var_soc)[1]/abs(current_data_range$var_soc - current_data$var_soc)[2:7])

abs(current_data_range$var_stat - current_data$var_stat)
mean(abs(current_data_range$var_stat - current_data$var_stat)[1]/abs(current_data_range$var_stat - current_data$var_stat)[2:6])


line_label = c('Microbial CUE', 'Carbon transfer fraction', 'Environmental modifer', 'Carbon input', 'Baseline decomposition', 'Vertical transport', 'Carbon input allocation')


S_sqrt = function(x){sign(x)*sqrt(abs(x))}
IS_sqrt = function(x){x^2*sign(x)}
S_sqrt_trans = function() trans_new("S_sqrt",S_sqrt,IS_sqrt)

p_sensitivity =
  ggplot(data = current_data) +
  geom_hline(yintercept = 0, size = 2, color = 'snow4', linetype = 'dashed', alpha = 1) +
  geom_vline(xintercept = 0, size = 2, color = 'snow4', linetype = 'dashed', alpha = 1) +
  geom_errorbar(aes(x = var_soc, ymin = lower_stat, ymax = upper_stat, color = component), size = 2, width = 0.7, stat = 'identity') +
  geom_errorbarh(aes(y = var_stat, xmin = lower_soc, xmax = upper_soc, color = component), size = 2, height = 0.02, stat = 'identity') +
  geom_point(aes(x = var_soc, y = var_stat, color = component, fill = component), shape = 15, size = 8) +
    geom_text(aes(x = var_soc, y = var_stat, label = line_label[order(component)], color = component), size = 10, vjust = c(0.5, 0.5, 3, 2, 0.5, -2, 0.5), hjust = c(1.2, -0.1, 0.2, 0.5, -0.1, 0.7, -0.1)) + # cue, I, K, V, Xi, NPP, transfer 
  scale_color_manual(name = '', labels = line_label, values = color_scheme) +
  scale_fill_manual(name = '', labels = line_label, values = color_scheme) +
  scale_x_continuous(limits = c(-20, 1200), breaks = c(0, 20, 90, 200, 400, 700, 1000), trans = 'S_sqrt') +
  scale_y_continuous(limits = c(NA, NA), breaks = -c(0, 0.01, 0.04, 0.09, 0.16, 0.25, 0.4), trans = 'S_sqrt') +
  theme_classic() +
  # change the legend properties
  # theme(panel.background = element_rect(fill = 'grey98'), plot.background = element_rect(fill = 'grey98'))+
  theme(legend.justification = 'right', legend.position = 'none', legend.background = element_rect(fill = NA), legend.text.align = 0) +
  theme(legend.text = element_text(size = 30), legend.title = element_text(size = 30))  +
  theme(legend.key = element_rect(color = NA, fill = NA), legend.key.size = unit(0.8, 'inch')) +
  # add title
  labs(x = 'Total absolute error of SOC estimates (Pg C)', y = expression(paste('Deviation of modeling efficiency'), sep = '')) +
  # modify the position of title
  # modify the font size
  theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), 'inch'), plot.background = element_rect(fill = 'transparent', color = NA)) +
  theme(axis.text = element_text(size = 30), axis.title = element_text(size = 35), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1), axis.ticks.length = unit(0.12, 'inch'))


#################################################################################
# sensitivity curve
#################################################################################
process_list =  c('A', 'K', 'Xi', 'V', 'I', 'E', 'NPP')
process_label = c('A_A', 'D_K', 'B_Xi', 'E_V', 'F_I', 'A_E', 'C_NPP')

process_description = c('Microbial CUE', 'Baseline Decomposition', 'Environmental Impacts', 'Vertical Transport', 'Input Allocation', 'Carbon transfer fraction', 'Carbon Input')
climate_list = c('A', 'B', 'C', 'D', 'E_all')

manage_proposal = seq(from = -0.2, to = 0.2, by = 0.04)*100
depth_name = c('C_30cm', 'B_100cm', 'A_200cm', 'D_full_depth')

marginal_change_total = readMat(paste(data_dir_output, 'world_simulation_analyses/marginal_change_summary_bootstrap_', model_name, '.mat', sep = ''))
marginal_change_total = marginal_change_total$marginal.change.total
####
marginal_change_total_mean = apply(marginal_change_total, c(1, 2, 3), mean, na.rm = TRUE)
marginal_change_total_upper = apply(marginal_change_total, c(1, 2, 3), quantile, probs = 0.025, na.rm = TRUE)
marginal_change_total_lower = apply(marginal_change_total, c(1, 2, 3), quantile, probs = 0.975, na.rm = TRUE)

col_list = c('stock', 'stock_upper', 'stock_lower', 'process', 'process_upper', 'process_lower', 'process_name')
marginal_change_summary = data.frame(array(NA, dim = c(length(process_list)*length(manage_proposal), length(col_list))))
colnames(marginal_change_summary) = col_list

counter = 1

for (iprocess in 1:length(process_list)) {
  
  for (imanage in 1:length(manage_proposal)) {
    marginal_change_summary[counter, 'process_name'] = process_label[iprocess]
    marginal_change_summary[counter, 'stock'] = marginal_change_total_mean[1, iprocess, imanage]
    marginal_change_summary[counter, 'stock_upper'] = marginal_change_total_upper[1, iprocess, imanage]
    marginal_change_summary[counter, 'stock_lower'] = marginal_change_total_lower[1, iprocess, imanage]
    
    marginal_change_summary[counter, 'process'] = marginal_change_total_mean[3, iprocess, imanage]
    marginal_change_summary[counter, 'process_upper'] = marginal_change_total_upper[3, iprocess, imanage]
    marginal_change_summary[counter, 'process_lower'] = marginal_change_total_lower[3, iprocess, imanage]
    
    counter = counter + 1
  }
}

#-----------------------------------------------------------------

stock_reference = marginal_change_summary$stock[6]
marginal_change_summary$stock = (marginal_change_summary$stock/stock_reference-1)*100
marginal_change_summary$stock_upper = (marginal_change_summary$stock_upper/stock_reference-1)*100
marginal_change_summary$stock_lower = (marginal_change_summary$stock_lower/stock_reference-1)*100

p_curve =
  ggplot(data = marginal_change_summary) +
  geom_ribbon(aes(x = process*100, ymin = stock_lower, ymax = stock_upper, fill = process_name), alpha = 0.15) +
  geom_line(aes(x = process*100, y = stock, color = process_name), alpha = 1, size = 2) +
  geom_hline(yintercept = 0, color = 'grey4', size = 2, linetype = 'dashed', alpha = 0.7) +
  # geom_hline(yintercept = marginal_change_summary$stock[6]*(1+0.004)**30, color = 'red', size = 2, linetype = 'dashed', alpha = 1) +
  geom_vline(xintercept = 0, color = 'grey4', size = 2, linetype = 'dashed', alpha = 0.7) +
  scale_x_continuous(trans = 'identity', limits = c(-10, 10), n.breaks = 7) +
  scale_y_continuous(trans = 'identity', limits = c(-30, 40), n.breaks = 7) +
  scale_color_manual(name = '', labels = line_label, values = color_scheme) +
  scale_fill_manual(name = '', labels = line_label, values = color_scheme) +
  # change the background to black and white
  theme_classic() +
  # theme(legend.position = 'None') +
  # theme(panel.background = element_rect(fill = 'grey98'), plot.background = element_rect(fill = 'grey98'))+
  theme(legend.justification = 'right', legend.position = 'right', legend.background = element_rect(fill = NA), legend.text.align = 0) +
  theme(legend.text = element_text(size = 30), legend.title = element_text(size = 30))  +
  theme(legend.key = element_rect(color = NA, fill = NA), legend.key.size = unit(0.8, 'inch')) +
  # add title
  labs(y = 'Changes of global SOC stock (%)', x = paste('Proportional changes of components (%)')) +
  # modify the position of title
  # modify the font sizea
  theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), 'inch'), plot.background = element_rect(fill = 'transparent', color = NA)) +
  theme(axis.text=element_text(size = 30), axis.title = element_text(size = 35), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1), axis.ticks.length = unit(0.12, 'inch'))


jpeg(paste('./Ensemble/revision4_main_fig4.jpeg', sep = ''), width = 26.5, height = 10, units = 'in', res = 300)
plot_grid(p_sensitivity, p_curve,
          nrow = 1, ncol = 2,
          rel_widths = c(1, 1.65),
          labels = c('a', 'b'),
          label_size = 50,
          label_x = 0.0, label_y = 1.02,
          label_fontfamily = 'Arial',
          label_fontface = 'bold'
)

dev.off()



pdf(paste('./Ensemble/revision4_main_fig4.pdf', sep = ''), width = 26.5, height = 10)

plot_grid(p_sensitivity, p_curve,
          nrow = 1, ncol = 2,
          rel_widths = c(1, 1.65),
          labels = c('a', 'b'),
          label_size = 50,
          label_x = 0.0, label_y = 1.02,
          label_fontface = 'bold'
)

dev.off()

###########################################
sensitivity_summary = array(NA, dim = c(6, 1))

iprocess = 5
for (iprocess in 1:length(process_list)) {
  
  regression_data_soc = (marginal_change_summary[marginal_change_summary$process_name== process_label[iprocess], 'stock']/100 + 1)*1
  regression_data_process = marginal_change_summary[marginal_change_summary$process_name== process_label[iprocess], 'process']
  
  
  lm_soc = lm(log(regression_data_soc/regression_data_soc[6]) ~ log(regression_data_process + 1))
  
  if (process_label[iprocess] == 'F_I') {
    
    lm_soc = lm(log((regression_data_soc/regression_data_soc[6])[8:11]) ~ log(regression_data_process[8:11] + 1))
  }
  
  # ggplot() +
  #   geom_line(aes(y = log(regression_data_soc/regression_data_soc[6]), x = log(regression_data_process + 1)), color = 'blue') +
  #   geom_line(aes(y = predict(lm_soc), x = log(regression_data_process + 1)), color = 'red')
  
  
  sensitivity_summary[iprocess] = exp((log((1 + 0.1)) - lm_soc$coefficients[1])/lm_soc$coefficients[2]) - 1
  
}

sensitivity_summary

