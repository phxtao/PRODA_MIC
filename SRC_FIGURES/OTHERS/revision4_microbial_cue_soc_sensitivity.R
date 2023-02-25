## Packages
library(R.matlab)
library(maps)
library(ggplot2)
library(R.matlab)
library(cowplot)
library(jcolors)
library(viridis)
library(cowplot)
library(ggnewscale)
# dev.off()
##
rm(list = ls())

setwd('/Users/phoenix/Google_Drive/R')

#############################################################################
# Data Path
#############################################################################
cue_soc_relation = readMat('/Users/phoenix/Google_Drive/Tsinghua_Luo/Projects/DATAHUB/ENSEMBLE/OUTPUT_DATA/mic_model_cue_soc_relationship/cue_soc_relation.mat')
cue_soc_relation = cue_soc_relation$cue.soc.relation

cue_soc_relation[cue_soc_relation < 0] = NA

tau_mic = c(0.1, 0.2, 0.5, 1, 1.5, 2, 5)
tau_enz = c(0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1)
mic_cue = seq(from = 0.1, to = 0.7, by = 0.1)



#############################different relationhsips
selected_mort = 2
selected_tau = 6

current_data_nort = c()
for (itau in 3:length(tau_enz)) {
  current_data_nort = rbind(current_data_nort, cbind(tau_mic[selected_mort], tau_enz[itau], mic_cue, cue_soc_relation[selected_mort, itau, ]/1000))
}

current_data_tau = c()
for (iallo in 3:length(tau_mic)) {
  current_data_tau = rbind(current_data_tau, cbind(tau_mic[iallo], tau_enz[selected_tau], mic_cue, cue_soc_relation[iallo, selected_tau, ]/1000))
}

current_data_nort = data.frame(current_data_nort)
colnames(current_data_nort) = c('allo', 'tau', 'cue', 'soc')

current_data_tau = data.frame(current_data_tau)
colnames(current_data_tau) = c('allo', 'tau', 'cue', 'soc')


jpeg(paste('./Ensemble/revision4_cue_soc_relation.jpeg', sep = ''), width = 13, height = 10, units = 'in', res = 300)

ggplot() + 
  geom_line(data = current_data_nort, aes(x = cue, y = soc, group = interaction(as.factor(allo), as.factor(tau)), color = as.factor(tau)), size = 2, linetype = 'solid') +
  geom_point(data = current_data_nort, aes(x = cue, y = soc, group = interaction(as.factor(allo), as.factor(tau)), color = as.factor(tau)), size = 6, shape = 16) +
  scale_color_manual(name = expression(paste(tau['ENZ,decay'], ' (', tau['MIC'], ' = 0.2)', sep = '')), values = rev(c('#08519c', '#3182bd', '#6baed6', '#9ecae1', '#c6dbef'))) +
  new_scale_color() +
  geom_line(data = current_data_tau, aes(x = cue, y = soc, group = interaction(as.factor(allo), as.factor(tau)), color = as.factor(allo)), size = 2, linetype = 'dashed') + 
  geom_point(data = current_data_tau, aes(x = cue, y = soc, group = interaction(as.factor(allo), as.factor(tau)), color = as.factor(allo)), size = 6, shape = 17) + 
  scale_color_manual(name = expression(paste(tau['MIC'], ' (', tau['ENZ,decay'], ' = 0.5)', sep = '')), values = rev(c('#a50f15', '#de2d26', '#fb6a4a', '#fc9272', '#fcbba1'))) +
  scale_x_continuous(trans = 'identity', limits = c(0.1, 0.7), breaks = mic_cue) +
  scale_y_continuous(trans = 'log10') + 
  # change the background to black and white
  theme_classic() +
  # add title
  labs(title = '', x = 'Carbon use efficiency', y = expression(paste('SOC stock (kg C m'^'-2', ')', sep = ''))) + 
  # modify the position of title
  theme(plot.title = element_text(hjust = 0.5, size = 50)) + 
  theme(legend.justification = c(0, 0), legend.position = 'right', legend.background = element_rect(fill = NA), legend.text.align = 0) +
  theme(legend.text = element_text(size = 25), legend.title = element_text(size = 25))  +
  theme(legend.key = element_rect(color = NA, fill = NA), legend.key.size = unit(0.5, 'inch')) +
  theme(legend.direction = 'vertical') + 
  # modify the margin
  # theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
  theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), 'inch')) +
  theme(axis.text=element_text(size = 30), axis.title = element_text(size = 35), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1), axis.ticks.length = unit(0.12, 'inch')) 

dev.off()
