% clc;
% clear;
%%
model_name = 'cesm2_clm5_mic_vr_v22';
nn_exp_name = 'exp_pc_cesm2_23';
time_domain = 'whole_time';

% server
data_dir_output = '/GFPS8p/cess11/taof/datahub/ensemble/output_data/';
data_dir_input = '/GFPS8p/cess11/taof/datahub/ensemble/input_data/';

% mac
% data_dir_output = '/Users/phoenix/Google_Drive/Tsinghua_Luo/Projects/DATAHUB/ENSEMBLE/OUTPUT_DATA/';
% data_dir_input = '/Users/phoenix/Google_Drive/Tsinghua_Luo/Projects/DATAHUB/ENSEMBLE/INPUT_DATA/';

%%
bootstrap_num = 200;

stat_matric = {'coeff_efficiency', 'correlation', 'coeff_determination', 'mse'};
stat_matric_index = 1;

component_name = {'H_NN', 'A_A', 'F_I', 'D_K', 'E_V', 'B_Xi', 'C_NPP', 'A_E', 'Ga_Homo'};

process_list =  {'A', 'K', 'Xi', 'V', 'I', 'E', 'NPP'};
process_label = {'A_A', 'D_K', 'B_Xi', 'E_V', 'F_I', 'A_E', 'C_NPP'};

process_description = {'CUE', 'Baseline Decomposition', 'Environmental Impacts', 'Vertical Transport', 'Input Allocation', 'Mortality Stabilization', 'Carbon Input'};
climate_list = {'A', 'B', 'C', 'D', 'E_all'};

manage_proposal = -0.1:0.02:0.1;
depth_name = {'C_30cm', 'B_100cm', 'A_200cm', 'D_full_depth'};

%% component control
control_summary = nan(length(depth_name)*length(component_name), 8);
colnames = {'var_soc', 'upper_soc', 'lower_soc', 'var_stat', 'upper_stat', 'lower_stat', 'depth', 'component'};

iquantile = 4;
ibootstrap = 1;
icomponent = 2;
counter = 1;
for iquantile = 1:4
  for icomponent = 1:length(component_name)
    control_matric_total_stock = nan(bootstrap_num, 1);
    control_matric_stat = nan(bootstrap_num, 1);
    
    for ibootstrap =1:bootstrap_num
      disp(['processing quantile ', num2str(iquantile), ' component ', num2str(icomponent), ' bootstrap  ', num2str(ibootstrap)]);
      % soc global stock
%       control_total_stock = load([data_dir_output, 'world_simulation_analyses/compontent_control_soc_total_stock_', model_name, '_', nn_exp_name, ...
%           '_cross_valid_0_', num2str(ibootstrap), '_control_test.mat']);
      control_total_stock = load([data_dir_output, 'world_simulation_analyses/compontent_control_soc_total_stock_', model_name, '_', nn_exp_name, ...
          '_bootstrap_', num2str(ibootstrap), '_control_test.mat']);
      
      control_total_stock = control_total_stock.var_data_middle;
      % record
      if icomponent == 1
        control_matric_total_stock(ibootstrap) = control_total_stock(iquantile, icomponent);
      else 
        control_matric_total_stock(ibootstrap) = control_total_stock(iquantile, icomponent);
      end
      
      
%       control_test = load([data_dir_output, 'world_simulation_analyses/compontent_control_spatial_variation_', model_name, '_', time_domain,...
%           '_', nn_exp_name, '_cross_valid_0_', num2str(ibootstrap), '.mat']);
      control_test = load([data_dir_output, 'world_simulation_analyses/compontent_control_spatial_variation_', model_name, '_', time_domain,...
          '_', nn_exp_name, '_bootstrap_', num2str(ibootstrap), '.mat']);
      
      control_test = control_test.var_data_middle; % [iquantile, icomponent, imatric]
      
      % record
	  if icomponent == 1
		control_matric_stat(ibootstrap) = control_test(icomponent, stat_matric_index);
      else
	  	control_matric_stat(ibootstrap) = control_test(icomponent, stat_matric_index) - control_test(1, stat_matric_index);
	  end
    end
    % summary
    control_summary(counter, 1) = median(control_matric_total_stock, 'omitnan');
    control_summary(counter, 2) = quantile(control_matric_total_stock, 0.975);
    control_summary(counter, 3) = quantile(control_matric_total_stock, 0.025);
    control_summary(counter, 4) = median(control_matric_stat, 'omitnan');
    control_summary(counter, 5) = quantile(control_matric_stat, 0.975);
    control_summary(counter, 6) = quantile(control_matric_stat, 0.025);
    
    
    control_summary(counter, 7) = iquantile;
    control_summary(counter, 8) = icomponent;
    
    counter = counter + 1;
    
  end
end

save([data_dir_output, 'world_simulation_analyses/compontent_control_summary_bootstrap_', model_name, '.mat'], 'control_summary');

%% sensitivity curve
grid_lon_lat = load([data_dir_output, 'world_simulation_analyses/marginal_sensitivity_proda_NPP_', model_name, '_', nn_exp_name, '_cross_valid_0_1.mat']);

grid_lon_lat = grid_lon_lat.var_data_middle{1};
grid_lon_lat = grid_lon_lat(:, 1:2);


% global soc stock
resolution = 0.5;
lat_seq = grid_lon_lat(: , 2);
% area of grid
radius = 6371008.8;
length_top = (2*pi*radius*cos(abs(lat_seq+resolution/2)/180*pi)/360)*resolution;
length_down = (2*pi*radius*cos(abs(lat_seq-resolution/2)/180*pi)/360)*resolution;
height = (pi*radius/180)*resolution;
lat_grid_area = (length_top + length_down)*height/2;

%-----------------------------------------------------------------
marginal_change_total = nan(3, length(process_list), length(manage_proposal), bootstrap_num);
uncertainty_process = nan(size(grid_lon_lat, 1), (length(process_list)+1), bootstrap_num);

counter = 1;
for iprocess = 1:length(process_list)
    bulk_process = nan(size(grid_lon_lat, 1), length(manage_proposal), bootstrap_num);
    soc_stock = nan(size(grid_lon_lat, 1), length(manage_proposal), bootstrap_num);
    
    for ibootstrap = 1:bootstrap_num
        disp(['processing process ', num2str(iprocess), ' bootstrap ', num2str(ibootstrap)]);
%         global_simu = load([data_dir_output, 'world_simulation_analyses/marginal_sensitivity_proda_', ...
%             process_list{iprocess}, '_', model_name, '_', nn_exp_name, '_cross_valid_0_', num2str(ibootstrap), '.mat']);
        global_simu = load([data_dir_output, 'world_simulation_analyses/marginal_sensitivity_proda_', ...
            process_list{iprocess}, '_', model_name, '_', nn_exp_name, '_bootstrap_', num2str(ibootstrap), '.mat']);
        global_simu = global_simu.var_data_middle;
        
        soc_stock(:, :, ibootstrap) = global_simu{2}; % soc stock
        bulk_process(:, :, ibootstrap) = global_simu{2+iprocess}; % corresponding managed process
       
		if iprocess == 1
			uncertainty_process(:, 1, ibootstrap) = soc_stock(:, 6, ibootstrap);
			uncertainty_process(:, (iprocess+1), ibootstrap) = bulk_process(:, 6, ibootstrap);								       
		elseif iprocess == 5
			uncertainty_process(:, (iprocess+1), ibootstrap) = bulk_process(:, 6, ibootstrap).^5;
		else
			uncertainty_process(:, (iprocess+1), ibootstrap) = bulk_process(:, 6, ibootstrap);
		end

        for imanage = 1:length(manage_proposal)
            invalid_loc = find(soc_stock(:, imanage, ibootstrap) < 0 ...
                | soc_stock(: , imanage, ibootstrap) > 1000000 ...
                | bulk_process(: , imanage, ibootstrap) < 0);
            soc_stock(invalid_loc, :, :) = nan;
            bulk_process(invalid_loc, :, :) = nan;
            
            
            % calculation of soc stock and tau
            % changes of soc stock
            marginal_change_total(1, iprocess, imanage, ibootstrap) =  ...
                sum(soc_stock(:, imanage, ibootstrap).*lat_grid_area, 'omitnan')/(10^15);
            % changes of bulk process
            if iprocess == 5
                marginal_change_total(3, iprocess, imanage, ibootstrap) = ...
                    median((bulk_process(: , imanage, ibootstrap).^5 - bulk_process(: , 6, ibootstrap).^5)./bulk_process(: , 6, ibootstrap).^5, 'omitnan');
                
            else
                marginal_change_total(3, iprocess, imanage, ibootstrap) = ...
                    median((bulk_process(: , imanage, ibootstrap) - bulk_process(: , 6, ibootstrap))./bulk_process(: , 6, ibootstrap), 'omitnan');
            end
        end
    end
end
uncertainty_process = std(uncertainty_process, [], 3, 'omitnan');
save([data_dir_output, 'world_simulation_analyses/uncertainty_process_bootstrap_', model_name, '.mat'], 'uncertainty_process');
save([data_dir_output, 'world_simulation_analyses/marginal_change_summary_bootstrap_', model_name, '.mat'], 'marginal_change_total');

