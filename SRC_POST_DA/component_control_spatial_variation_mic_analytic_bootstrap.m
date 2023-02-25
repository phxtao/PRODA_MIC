% clear;
% clc;
% start_id = 1;
% end_id = 1;
%%
% parallel setting
delete(gcp('nocreate'));
parpool(10);

cesm2_case_name = 'sasu_f05_g16_checked_step4';
start_year = 661;
end_year = 680;

model_name = 'cesm2_clm5_mic_vr_v22'; % 'cesm2_clm5_cen_vr_v2';

time_domain = 'whole_time'; % 'whole_time', 'before_1985', 'after_1985', 'random_half', 'random_half_2'


parfor ibootstrap = start_id:end_id
% for ibootstrap = start_id:10:end_id
    disp([datestr(now), ' Processing ibootstrap ', num2str(ibootstrap)]);
    
	% exp_name = ['exp_pc_cesm2_24_cross_valid_0_', num2str(ibootstrap)];
    exp_name = ['exp_pc_cesm2_24_bootstrap_', num2str(ibootstrap)];
    
    
    %% paths
    % mac
    % data_path = '/Users/phoenix/Google_Drive/Tsinghua_Luo/Projects/DATAHUB/ENSEMBLE/';
    % cd(['/Users/phoenix/Google_Drive/Tsinghua_Luo/Projects/ENSEMBLE/SRC_DA/src_', model_name, '/']);
    
    % server
    data_path = '/GFPS8p/cess11/taof/datahub/ensemble/';
    cd(['/GFPS8p/cess11/taof/ensemble/src_da/src_', model_name, '/']);
    
    %% set vertical soil pools
    month_num = 12;
    soil_cpool_num = 7;
    soil_decom_num = 20;
    
    %% load wosis data
    env_info = load([data_path, 'input_data/wosis_2019_snap_shot/wosis_2019_snapshot_hugelius_mishra_env_info.mat'], 'EnvInfo'); % NPP in this file will be used
    env_info = env_info.EnvInfo;
    % layer_info: "profile_id, date, upper_depth, lower_depth, node_depth, soc_layer_weight, soc_stock, bulk_denstiy, is_pedo"
    wosis_profile_info = ncread([data_path, 'input_data/wosis_2019_snap_shot/soc_profile_wosis_2019_snapshot_hugelius_mishra.nc'], 'soc_profile_info'); % wosis profile info
    wosis_soc_info = ncread([data_path, 'input_data/wosis_2019_snap_shot/soc_data_integrate_wosis_2019_snapshot_hugelius_mishra.nc'], 'data_soc_integrate'); % wosis SOC info
    
    %% soil depth information
    % width between two interfaces
    dz = [2.000000000000000E-002, 4.000000000000000E-002, 6.000000000000000E-002, ...
        8.000000000000000E-002, 0.120000000000000, 0.160000000000000, ...
        0.200000000000000, 0.240000000000000, 0.280000000000000, ...
        0.320000000000000, 0.360000000000000, 0.400000000000000, ...
        0.440000000000000, 0.540000000000000, 0.640000000000000, ...
        0.740000000000000, 0.840000000000000, 0.940000000000000, ...
        1.04000000000000, 1.14000000000000, 2.39000000000000, ...
        4.67553390593274, 7.63519052838329, 11.1400000000000, ...
        15.1154248593737]';
    
    % depth of the interface
    zisoi = [2.000000000000000E-002, 6.000000000000000E-002, ...
        0.120000000000000, 0.200000000000000, 0.320000000000000, ...
        0.480000000000000, 0.680000000000000, 0.920000000000000, ...
        1.20000000000000, 1.52000000000000, 1.88000000000000, ...
        2.28000000000000, 2.72000000000000, 3.26000000000000, ...
        3.90000000000000, 4.64000000000000, 5.48000000000000, ...
        6.42000000000000, 7.46000000000000, 8.60000000000000, ...
        10.9900000000000, 15.6655339059327, 23.3007244343160, ...
        34.4407244343160, 49.5561492936897]';
    
    % depth of the node
    zsoi = [1.000000000000000E-002, 4.000000000000000E-002, 9.000000000000000E-002, ...
        0.160000000000000, 0.260000000000000, 0.400000000000000, ...
        0.580000000000000, 0.800000000000000, 1.06000000000000, ...
        1.36000000000000, 1.70000000000000, 2.08000000000000, ...
        2.50000000000000, 2.99000000000000, 3.58000000000000, ...
        4.27000000000000, 5.06000000000000, 5.95000000000000, ...
        6.94000000000000, 8.03000000000000, 9.79500000000000, ...
        13.3277669529664, 19.4831291701244, 28.8707244343160, ...
        41.9984368640029]';
    
    % depth between two node
    dz_node = zsoi - [0; zsoi(1:end-1)];
    
    %% input from cesm2
    % cesm2 resolution
    cesm2_resolution_lat = 180/384;
    cesm2_resolution_lon = 360/576;
    lon_grid = (-180 + cesm2_resolution_lon/2 : cesm2_resolution_lon : 180 - cesm2_resolution_lon/2)';
    lat_grid = (90 - cesm2_resolution_lat/2 : -cesm2_resolution_lat : -90 + cesm2_resolution_lat/2)';
    
    % nbedrock
    cesm2_simu_nbedrock = load([data_path, 'input_data/cesm2_simu/spinup_ss/',...
        cesm2_case_name, '_cesm2_ss_4da_', num2str(start_year), '_', num2str(end_year),...
        '_nbedrock.mat'], 'var_record_monthly_mean');
    cesm2_simu_nbedrock = cesm2_simu_nbedrock.var_record_monthly_mean;
    % ALTMAX
    cesm2_simu_altmax = load([data_path, 'input_data/cesm2_simu/spinup_ss/',...
        cesm2_case_name, '_cesm2_ss_4da_', num2str(start_year), '_', num2str(end_year),...
        '_ALTMAX.mat'], 'var_record_monthly_mean');
    cesm2_simu_altmax = cesm2_simu_altmax.var_record_monthly_mean;
    % ALTMAX_LASTYEAR
    cesm2_simu_altmax_last_year = load([data_path, 'input_data/cesm2_simu/spinup_ss/',...
        cesm2_case_name, '_cesm2_ss_4da_', num2str(start_year), '_', num2str(end_year),...
        '_ALTMAX_LASTYEAR.mat'], 'var_record_monthly_mean');
    cesm2_simu_altmax_last_year = cesm2_simu_altmax_last_year.var_record_monthly_mean;
    % CELLSAND
    cesm2_simu_cellsand = load([data_path, 'input_data/cesm2_simu/spinup_ss/',...
        cesm2_case_name, '_cesm2_ss_4da_', num2str(start_year), '_', num2str(end_year),...
        '_CELLSAND.mat'], 'var_record_monthly_mean');
    cesm2_simu_cellsand = cesm2_simu_cellsand.var_record_monthly_mean;
    % NPP
    cesm2_simu_npp = load([data_path, 'input_data/cesm2_simu/spinup_ss/',...
        cesm2_case_name, '_cesm2_ss_4da_', num2str(start_year), '_', num2str(end_year),...
        '_NPP.mat'], 'var_record_monthly_mean');
    cesm2_simu_npp = cesm2_simu_npp.var_record_monthly_mean;
    % SOILPSI
    cesm2_simu_soil_water_potnetial = load([data_path, 'input_data/cesm2_simu/spinup_ss/',...
        cesm2_case_name, '_cesm2_ss_4da_', num2str(start_year), '_', num2str(end_year),...
        '_SOILPSI.mat'], 'var_record_monthly_mean');
    cesm2_simu_soil_water_potnetial = cesm2_simu_soil_water_potnetial.var_record_monthly_mean;
    % TSOI
    cesm2_simu_soil_temperature = load([data_path, 'input_data/cesm2_simu/spinup_ss/',...
        cesm2_case_name, '_cesm2_ss_4da_', num2str(start_year), '_', num2str(end_year),...
        '_TSOI.mat'], 'var_record_monthly_mean');
    cesm2_simu_soil_temperature = cesm2_simu_soil_temperature.var_record_monthly_mean;
    % W_SCALAR
    cesm2_simu_w_scalar = load([data_path, 'input_data/cesm2_simu/spinup_ss/',...
        cesm2_case_name, '_cesm2_ss_4da_', num2str(start_year), '_', num2str(end_year),...
        '_W_SCALAR.mat'], 'var_record_monthly_mean');
    cesm2_simu_w_scalar = cesm2_simu_w_scalar.var_record_monthly_mean;
    % T_SCALAR
    cesm2_simu_t_scalar = load([data_path, 'input_data/cesm2_simu/spinup_ss/',...
        cesm2_case_name, '_cesm2_ss_4da_', num2str(start_year), '_', num2str(end_year),...
        '_T_SCALAR.mat'], 'var_record_monthly_mean');
    cesm2_simu_t_scalar = cesm2_simu_t_scalar.var_record_monthly_mean;
    % O_SCALAR
    cesm2_simu_o_scalar = load([data_path, 'input_data/cesm2_simu/spinup_ss/',...
        cesm2_case_name, '_cesm2_ss_4da_', num2str(start_year), '_', num2str(end_year),...
        '_O_SCALAR.mat'], 'var_record_monthly_mean');
    cesm2_simu_o_scalar = cesm2_simu_o_scalar.var_record_monthly_mean;
    % FPI_vr
    cesm2_simu_n_scalar = load([data_path, 'input_data/cesm2_simu/spinup_ss/',...
        cesm2_case_name, '_cesm2_ss_4da_', num2str(start_year), '_', num2str(end_year),...
        '_FPI_vr.mat'], 'var_record_monthly_mean');
    cesm2_simu_n_scalar = cesm2_simu_n_scalar.var_record_monthly_mean;
    % LITR1_INPUT_ACC_VECTOR
    cesm2_simu_input_vector_litter1 = load([data_path, 'input_data/cesm2_simu/spinup_ss/',...
        cesm2_case_name, '_cesm2_ss_4da_', num2str(start_year), '_', num2str(end_year),...
        '_LITR1_INPUT_ACC_VECTOR.mat'], 'var_record_monthly_mean');
    cesm2_simu_input_vector_litter1 = cesm2_simu_input_vector_litter1.var_record_monthly_mean;
    % LITR2_INPUT_ACC_VECTOR
    cesm2_simu_input_vector_litter2 = load([data_path, 'input_data/cesm2_simu/spinup_ss/',...
        cesm2_case_name, '_cesm2_ss_4da_', num2str(start_year), '_', num2str(end_year),...
        '_LITR2_INPUT_ACC_VECTOR.mat'], 'var_record_monthly_mean');
    cesm2_simu_input_vector_litter2 = cesm2_simu_input_vector_litter2.var_record_monthly_mean;
    % LITR3_INPUT_ACC_VECTOR
    cesm2_simu_input_vector_litter3 = load([data_path, 'input_data/cesm2_simu/spinup_ss/',...
        cesm2_case_name, '_cesm2_ss_4da_', num2str(start_year), '_', num2str(end_year),...
        '_LITR3_INPUT_ACC_VECTOR.mat'], 'var_record_monthly_mean');
    cesm2_simu_input_vector_litter3 = cesm2_simu_input_vector_litter3.var_record_monthly_mean;
    % CWD_INPUT_ACC_VECTOR
    cesm2_simu_input_vector_cwd = load([data_path, 'input_data/cesm2_simu/spinup_ss/',...
        cesm2_case_name, '_cesm2_ss_4da_', num2str(start_year), '_', num2str(end_year),...
        '_CWD_INPUT_ACC_VECTOR.mat'], 'var_record_monthly_mean');
    cesm2_simu_input_vector_cwd = cesm2_simu_input_vector_cwd.var_record_monthly_mean;
    % TOTSOMC
    cesm2_simu_soc_stock = load([data_path, 'input_data/cesm2_simu/spinup_ss/',...
        cesm2_case_name, '_cesm2_ss_4da_', num2str(start_year), '_', num2str(end_year),...
        '_TOTSOMC.mat'], 'var_record_monthly_mean');
    cesm2_simu_soc_stock = cesm2_simu_soc_stock.var_record_monthly_mean;
    
    
    for ilayer = 1:soil_decom_num
        cesm2_simu_input_vector_litter1(:, :, ilayer, :) = cesm2_simu_input_vector_litter1(:, :, ilayer, :).*dz(ilayer);
        cesm2_simu_input_vector_litter2(:, :, ilayer, :) = cesm2_simu_input_vector_litter2(:, :, ilayer, :).*dz(ilayer);
        cesm2_simu_input_vector_litter3(:, :, ilayer, :) = cesm2_simu_input_vector_litter3(:, :, ilayer, :).*dz(ilayer);
        cesm2_simu_input_vector_cwd(:, :, ilayer, :) = cesm2_simu_input_vector_cwd(:, :, ilayer, :).*dz(ilayer);
    end
    
    cesm2_simu_input_sum_litter1 = reshape(sum(cesm2_simu_input_vector_litter1, 3), [384, 576, 12]);
    cesm2_simu_input_sum_litter2 = reshape(sum(cesm2_simu_input_vector_litter2, 3), [384, 576, 12]);
    cesm2_simu_input_sum_litter3 = reshape(sum(cesm2_simu_input_vector_litter3, 3), [384, 576, 12]);
    cesm2_simu_input_sum_cwd = reshape(sum(cesm2_simu_input_vector_cwd, 3), [384, 576, 12]);
    
    
    % clearvars cesm2_simu_input_vector_litter1 cesm2_simu_input_vector_litter2 cesm2_simu_input_vector_litter3 cesm2_simu_input_vector_cwd
    
    %% Difference between NN and MCMC results
    para_mean = load([data_path, 'output_data/mcmc_summary_', model_name, '/', model_name '_para_mean.mat']);
    para_mean = para_mean.para_mean;
    stat_r2 = load([data_path, 'output_data/mcmc_summary_', model_name, '/', model_name, '_stat_r2.mat']);
    stat_r2 = stat_r2.stat_r2;
    stat_r2 = max(stat_r2, [], 2);
    
    stat_gr = load([data_path, 'output_data/mcmc_summary_', model_name, '/', model_name, '_para_gr.mat']);
    stat_gr = stat_gr.para_gr;
    stat_gr = mean(stat_gr, 2);
    
    %% load results from neural networking
    nn_predict = csvread([data_path, 'output_data/neural_networking/nn_para_result_', model_name, '_', time_domain, '_', exp_name, '.csv']);
    nn_site_loc = csvread([data_path, 'output_data/neural_networking/nn_site_loc_', model_name, '_', time_domain, '_', exp_name, '.csv']);
    %     nn_site_loc = find(stat_r2 > 0 & stat_gr < 1.05);
    %     nn_predict = para_mean(nn_site_loc, :);
    
    para_global_constant = mean(para_mean, 1, 'omitnan')';
    %     para_global_constant = mean(para_mean(randsample(length(nn_site_loc), 100, 'true'), :), 1, 'omitnan')';
    
    %% Individual Simulation
    indi_map_project_mod = nan(length(nn_site_loc), 150, 9);
    indi_map_project_obs = nan(length(nn_site_loc), 150, 9);
    indi_map_project_depth = nan(length(nn_site_loc), 150, 9);
    
    site_by_site_mod = nan(length(nn_site_loc), soil_decom_num);
    
    for iprofile_hat = 1:length(nn_site_loc)
        
        iprofile = nn_site_loc(iprofile_hat);
        warning('off')
        
        % disp([datestr(now), ' Processing profile ', num2str(iprofile_hat)]);
        
        profile_id = wosis_profile_info(iprofile, 1);
        % find currently using profile
        loc_profile = find(wosis_soc_info(:, 1) == profile_id);
        % info of the node depth of profile, and change unit from cm to m
        wosis_layer_depth = wosis_soc_info(loc_profile, 5)/100;
        % observed C info (gC/m3)
        wosis_layer_obs = wosis_soc_info(loc_profile, 7);
        % specify the number of layers in studied profiles
        num_layers = length(wosis_layer_obs);
        
        % find the lon and lat info of soil profile
        lon_profile = wosis_profile_info(iprofile, 4);
        lat_profile = wosis_profile_info(iprofile, 5);
        lat_loc = find(abs(lat_profile - lat_grid) == min(abs(lat_profile - lat_grid)));
        lon_loc = find(abs(lon_profile - lon_grid) == min(abs(lon_profile - lon_grid)));
        
        if length(lon_loc) > 1
            lon_loc = lon_loc(1);
        end
        
        if length(lat_loc) > 1
            lat_loc = lat_loc(1);
        end
        
        % input vector
        
        input_vector_cwd = reshape(cesm2_simu_input_sum_cwd(lat_loc, lon_loc, :), [1, month_num]);
        input_vector_litter1 = reshape(cesm2_simu_input_sum_litter1(lat_loc, lon_loc, :), [1, month_num]);
        input_vector_litter2 = reshape(cesm2_simu_input_sum_litter2(lat_loc, lon_loc, :), [1, month_num]);
        input_vector_litter3 = reshape(cesm2_simu_input_sum_litter3(lat_loc, lon_loc, :), [1, month_num]);
        
        % no input information
        if isnan(input_vector_cwd(1)) == 1 || isnan(input_vector_litter1(1)) == 1 ...
                || isnan(input_vector_litter2(1)) == 1 || isnan(input_vector_litter3(1)) == 1
            disp(['No input info ', ' Profile ', num2str(iprofile)]);
            continue
        end
        
        % npp from CESM2 simulatoin
        npp_mean = reshape(cesm2_simu_npp(lat_loc, lon_loc, :), [month_num, 1]);
        % NPP info of studied profile from soc_modIS NPP mean (year 2001-2016)
        % modis_npp_mean = env_info(iprofile, 15);
        
        % altmax current and last year
        altmax_current_profile = reshape(cesm2_simu_altmax(lat_loc, lon_loc, :), [month_num, 1]);
        altmax_lastyear_profile = reshape(cesm2_simu_altmax_last_year(lat_loc, lon_loc, :), [month_num, 1]);
        
        % nbedrock
        nbedrock = reshape(cesm2_simu_nbedrock(lat_loc, lon_loc, :), [month_num, 1]);
        
        % oxygen scalar
        xio = reshape(cesm2_simu_o_scalar(lat_loc, lon_loc, 1:soil_decom_num, :), [soil_decom_num, month_num, 1]);
        % nitrogen
        xin = reshape(cesm2_simu_n_scalar(lat_loc, lon_loc, 1:soil_decom_num, :), [soil_decom_num, month_num, 1]);
        
        % sand content
        sand = reshape(cesm2_simu_cellsand(lat_loc, lon_loc, 1:soil_decom_num, :), [soil_decom_num, month_num, 1]);
        sand_vector = sand;
        
        % soil temperature and water potential
        soil_temp_profile = reshape(cesm2_simu_soil_temperature(lat_loc, lon_loc, 1:soil_decom_num, :), [soil_decom_num, month_num, 1]);
        % soil_water_profile = reshape(cesm2_simu_soil_water_potnetial(lat_loc, lon_loc, 1:soil_decom_num, :), [soil_decom_num, month_num, 1]);
        % soil_temp_profile = reshape(cesm2_simu_t_scalar(lat_loc, lon_loc, 1:soil_decom_num, :), [soil_decom_num, month_num, 1]);
        soil_water_profile = reshape(cesm2_simu_w_scalar(lat_loc, lon_loc, 1:soil_decom_num, :), [soil_decom_num, month_num, 1]);
        
        %% data cleansing
        % deal with zero value in mcmc
        % if zero, cost function will always be inf
        if length(find(wosis_layer_obs == 0)) == length(wosis_layer_obs)
            disp(['All zero value ', ' Profile ', num2str(iprofile)]);
            continue
        elseif isempty(find(wosis_layer_obs == 0, 1)) == 0
            wosis_layer_depth = wosis_layer_depth(wosis_layer_obs > 0);
            wosis_layer_obs = wosis_layer_obs(wosis_layer_obs > 0);
        end
        
        % deal with nan values in obs or depth
        layer_valid_loc = find(isnan(wosis_layer_depth) == 0 & isnan(wosis_layer_obs) == 0);
        
        if isempty(layer_valid_loc) == 0
            wosis_layer_depth = wosis_layer_depth(layer_valid_loc);
            wosis_layer_obs = wosis_layer_obs(layer_valid_loc);
        else
            disp(['No valid obsrvations in ', ' Profile ', num2str(iprofile)]);
            continue
        end
        
        % eliminate profiles having single layer
        if length(wosis_layer_obs) == 1
            disp(['Single layer ', ' Profile ', num2str(iprofile)]);
            continue
        end
        
        %% Parameterization
        para_name = {'bio', 'cryo', 'q10', 'w_scaling', ...
            'tau4cwd', 'tau4l1', 'tau4l2', 'tau4s1', 'tau4s2_death', 'tau4s2_enz', 'tau4s3', 'tau4s4', ...
            'mm_const_assim', 'mm_const_decom', ...
            'fcwdl2', ...
            'pl1s1', 'pl2s1', 'pl3s4', 'l1_cue', 'l2l3_cue', ...
            'mic_cue', 'pdeath2soc', ...
            'beta'};
        
        para_A_loc = [8, 13, 19:21];
        para_K_loc = [5:7, 9:12, 14];
        % para_U_loc = [8, 13];
        para_Xi_loc = [3, 4];
        para_I_loc = 23;
        para_V_loc = [1, 2];
        para_E_loc = [15:18, 22];
                        
        %% Fully Varying (NN)
        control_loc = 1;
        % parameters values
        para = nn_predict(nn_site_loc == iprofile, :);
        % clear warning info
        lastwarn('');
        soc_mod = ...
            fun_var_decom(para, nbedrock, sand_vector, npp_mean, ...
            input_vector_cwd, input_vector_litter1, input_vector_litter2, input_vector_litter3, ...
            altmax_current_profile, altmax_lastyear_profile, soil_temp_profile, soil_water_profile, xio, xin);
        soc_mod = soc_mod(:, 5);
        
        [~, msgid] = lastwarn;
        if strcmp(msgid, 'MATLAB:singularMatrix') || strcmp(msgid, 'MATLAB:nearlySingularMatrix')
            disp(['Assilimation_failed_in_Profile_', num2str(iprofile)]);
            % msgid = [];
            continue
        end
        
        % profiles
        if length(wosis_layer_obs) > 1 && isnan(soc_mod(1)) == 0
            optimize_profile_soc = interp1(zsoi(1:soil_decom_num, 1), soc_mod, wosis_layer_depth, 'pchip');
            % scatter(current_obsC(:,1), current_obsC(:,2));
            site_by_site_mod(iprofile_hat, :) = soc_mod;
            indi_map_project_mod(iprofile_hat, 1:length(wosis_layer_obs), control_loc) = optimize_profile_soc;
            indi_map_project_obs(iprofile_hat, 1:length(wosis_layer_obs), control_loc) = wosis_layer_obs;
            indi_map_project_depth(iprofile_hat, 1:length(wosis_layer_obs), control_loc) = wosis_layer_depth;
        end
        
        
        %% component A
        control_loc = 2;
        % parameters values
        para = nn_predict(iprofile_hat, :);
        para(para_A_loc) = para_global_constant(para_A_loc);
        
        % clear warning info
        lastwarn('');
        
        soc_mod = ...
            fun_var_decom(para, nbedrock, sand_vector, npp_mean, ...
            input_vector_cwd, input_vector_litter1, input_vector_litter2, input_vector_litter3, ...
            altmax_current_profile, altmax_lastyear_profile, soil_temp_profile, soil_water_profile, xio, xin);
        soc_mod = soc_mod(:, 5);
        
        
        [~, msgid] = lastwarn;
        if strcmp(msgid, 'MATLAB:singularMatrix') || strcmp(msgid, 'MATLAB:nearlySingularMatrix')
            disp(['Assilimation_failed_in_Profile_', num2str(iprofile)]);
            % msgid = [];
            continue
        end
        
        % profiles
        if length(wosis_layer_obs) > 1 && isnan(soc_mod(1)) == 0
            optimize_profile_soc = interp1(zsoi(1:soil_decom_num, 1), soc_mod, wosis_layer_depth, 'pchip');
            % scatter(current_obsC(:,1), current_obsC(:,2));
            site_by_site_mod(iprofile_hat, :) = soc_mod;
            indi_map_project_mod(iprofile_hat, 1:length(wosis_layer_obs), control_loc) = optimize_profile_soc;
            indi_map_project_obs(iprofile_hat, 1:length(wosis_layer_obs), control_loc) = wosis_layer_obs;
            indi_map_project_depth(iprofile_hat, 1:length(wosis_layer_obs), control_loc) = wosis_layer_depth;
            
        end
        
        %% component I
        control_loc = 3;
        % parameters values
        para = nn_predict(iprofile_hat, :);
        para(para_I_loc) = para_global_constant(para_I_loc);
        
        % clear warning info
        lastwarn('');
        soc_mod = ...
            fun_var_decom(para, nbedrock, sand_vector, npp_mean, ...
            input_vector_cwd, input_vector_litter1, input_vector_litter2, input_vector_litter3, ...
            altmax_current_profile, altmax_lastyear_profile, soil_temp_profile, soil_water_profile, xio, xin);
        soc_mod = soc_mod(:, 5);
        
        [~, msgid] = lastwarn;
        if strcmp(msgid, 'MATLAB:singularMatrix') || strcmp(msgid, 'MATLAB:nearlySingularMatrix')
            disp(['Assilimation_failed_in_Profile_', num2str(iprofile)]);
            % msgid = [];
            continue
        end
        
        % profiles
        if length(wosis_layer_obs) > 1 && isnan(soc_mod(1)) == 0
            optimize_profile_soc = interp1(zsoi(1:soil_decom_num, 1), soc_mod, wosis_layer_depth, 'pchip');
            % scatter(current_obsC(:,1), current_obsC(:,2));
            site_by_site_mod(iprofile_hat, :) = soc_mod;
            indi_map_project_mod(iprofile_hat, 1:length(wosis_layer_obs), control_loc) = optimize_profile_soc;
            indi_map_project_obs(iprofile_hat, 1:length(wosis_layer_obs), control_loc) = wosis_layer_obs;
            indi_map_project_depth(iprofile_hat, 1:length(wosis_layer_obs), control_loc) = wosis_layer_depth;
            
        end
        
        %% component K
        control_loc = 4;
        % parameters values
        para = nn_predict(iprofile_hat, :);
        para(para_K_loc) = para_global_constant(para_K_loc);
        
        % clear warning info
        lastwarn('');
        soc_mod = ...
            fun_var_decom(para, nbedrock, sand_vector, npp_mean, ...
            input_vector_cwd, input_vector_litter1, input_vector_litter2, input_vector_litter3, ...
            altmax_current_profile, altmax_lastyear_profile, soil_temp_profile, soil_water_profile, xio, xin);
        soc_mod = soc_mod(:, 5);
        
        [~, msgid] = lastwarn;
        if strcmp(msgid, 'MATLAB:singularMatrix') || strcmp(msgid, 'MATLAB:nearlySingularMatrix')
            disp(['Assilimation_failed_in_Profile_', num2str(iprofile)]);
            % msgid = [];
            continue
        end
        
        % profiles
        if length(wosis_layer_obs) > 1 && isnan(soc_mod(1)) == 0
            optimize_profile_soc = interp1(zsoi(1:soil_decom_num, 1), soc_mod, wosis_layer_depth, 'pchip');
            % scatter(current_obsC(:,1), current_obsC(:,2));
            site_by_site_mod(iprofile_hat, :) = soc_mod;
            indi_map_project_mod(iprofile_hat, 1:length(wosis_layer_obs), control_loc) = optimize_profile_soc;
            indi_map_project_obs(iprofile_hat, 1:length(wosis_layer_obs), control_loc) = wosis_layer_obs;
            indi_map_project_depth(iprofile_hat, 1:length(wosis_layer_obs), control_loc) = wosis_layer_depth;
            
        end
        
        
%         %% component U
%         control_loc = 5;
%         % parameters values
%         para = nn_predict(iprofile_hat, :);
%         para(para_U_loc) = para_global_constant(para_U_loc);
%         
%         % clear warning info
%         lastwarn('');
%         soc_mod = ...
%             fun_var_decom(para, nbedrock, sand_vector, npp_mean, ...
%             input_vector_cwd, input_vector_litter1, input_vector_litter2, input_vector_litter3, ...
%             altmax_current_profile, altmax_lastyear_profile, soil_temp_profile, soil_water_profile, xio, xin);
%         soc_mod = soc_mod(:, 5);
%         
%         [~, msgid] = lastwarn;
%         if strcmp(msgid, 'MATLAB:singularMatrix') || strcmp(msgid, 'MATLAB:nearlySingularMatrix')
%             disp(['Assilimation_failed_in_Profile_', num2str(iprofile)]);
%             % msgid = [];
%             continue
%         end
%         
%         % profiles
%         if length(wosis_layer_obs) > 1 && isnan(soc_mod(1)) == 0
%             optimize_profile_soc = interp1(zsoi(1:soil_decom_num, 1), soc_mod, wosis_layer_depth, 'pchip');
%             % scatter(current_obsC(:,1), current_obsC(:,2));
%             site_by_site_mod(iprofile_hat, :) = soc_mod;
%             indi_map_project_mod(iprofile_hat, 1:length(wosis_layer_obs), control_loc) = optimize_profile_soc;
%             indi_map_project_obs(iprofile_hat, 1:length(wosis_layer_obs), control_loc) = wosis_layer_obs;
%             indi_map_project_depth(iprofile_hat, 1:length(wosis_layer_obs), control_loc) = wosis_layer_depth;
%             
%         end
        

        %% component V
        control_loc = 5;
        % parameters values
        para = nn_predict(iprofile_hat, :);
        para(para_V_loc) = para_global_constant(para_V_loc);
        
        % clear warning info
        lastwarn('');
        soc_mod = ...
            fun_var_decom(para, nbedrock, sand_vector, npp_mean, ...
            input_vector_cwd, input_vector_litter1, input_vector_litter2, input_vector_litter3, ...
            altmax_current_profile, altmax_lastyear_profile, soil_temp_profile, soil_water_profile, xio, xin);
        soc_mod = soc_mod(:, 5);
        
        [~, msgid] = lastwarn;
        if strcmp(msgid, 'MATLAB:singularMatrix') || strcmp(msgid, 'MATLAB:nearlySingularMatrix')
            disp(['Assilimation_failed_in_Profile_', num2str(iprofile)]);
            % msgid = [];
            continue
        end
        
        % profiles
        if length(wosis_layer_obs) > 1 && isnan(soc_mod(1)) == 0
            optimize_profile_soc = interp1(zsoi(1:soil_decom_num, 1), soc_mod, wosis_layer_depth, 'pchip');
            % scatter(current_obsC(:,1), current_obsC(:,2));
            site_by_site_mod(iprofile_hat, :) = soc_mod;
            indi_map_project_mod(iprofile_hat, 1:length(wosis_layer_obs), control_loc) = optimize_profile_soc;
            indi_map_project_obs(iprofile_hat, 1:length(wosis_layer_obs), control_loc) = wosis_layer_obs;
            indi_map_project_depth(iprofile_hat, 1:length(wosis_layer_obs), control_loc) = wosis_layer_depth;
            
        end
        
        %% component Xi
        control_loc = 6;
        % parameters values
        para = nn_predict(iprofile_hat, :);
        para(para_Xi_loc) = para_global_constant(para_Xi_loc);
        
        % clear warning info
        lastwarn('');
        soc_mod = ...
            fun_var_decom(para, nbedrock, sand_vector, npp_mean, ...
            input_vector_cwd, input_vector_litter1, input_vector_litter2, input_vector_litter3, ...
            altmax_current_profile, altmax_lastyear_profile, soil_temp_profile, soil_water_profile, xio, xin);
        soc_mod = soc_mod(:, 5);
        
        [~, msgid] = lastwarn;
        if strcmp(msgid, 'MATLAB:singularMatrix') || strcmp(msgid, 'MATLAB:nearlySingularMatrix')
            disp(['Assilimation_failed_in_Profile_', num2str(iprofile)]);
            % msgid = [];
            continue
        end
        
        % profiles
        if length(wosis_layer_obs) > 1 && isnan(soc_mod(1)) == 0
            optimize_profile_soc = interp1(zsoi(1:soil_decom_num, 1), soc_mod, wosis_layer_depth, 'pchip');
            % scatter(current_obsC(:,1), current_obsC(:,2));
            site_by_site_mod(iprofile_hat, :) = soc_mod;
            indi_map_project_mod(iprofile_hat, 1:length(wosis_layer_obs), control_loc) = optimize_profile_soc;
            indi_map_project_obs(iprofile_hat, 1:length(wosis_layer_obs), control_loc) = wosis_layer_obs;
            indi_map_project_depth(iprofile_hat, 1:length(wosis_layer_obs), control_loc) = wosis_layer_depth;
            
        end
        
        %% component NPP
        control_loc = 7;
        % parameters values
        para = nn_predict(iprofile_hat, :);
        
        days_in_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]';
        
        while true
            npp_scale_factor = 1 + randi([-1, 1])*2*(env_info(iprofile, 76)/env_info(iprofile, 75));
            if npp_scale_factor ~= 1
                break
            end
        end
        
        % npp_scale_factor = min(max(env_info(iprofile, 15)/env_info(iprofile, 75), 0.5), 1.5);
        % clear warning info
        lastwarn('');
        soc_mod = ...
            fun_var_decom(para, nbedrock, sand_vector, npp_mean*npp_scale_factor, ...
            input_vector_cwd*npp_scale_factor, input_vector_litter1*npp_scale_factor, input_vector_litter2*npp_scale_factor, input_vector_litter3*npp_scale_factor, ...
            altmax_current_profile, altmax_lastyear_profile, soil_temp_profile, soil_water_profile, xio, xin);
        soc_mod = soc_mod(:, 5);
        
        [~, msgid] = lastwarn;
        if strcmp(msgid, 'MATLAB:singularMatrix') || strcmp(msgid, 'MATLAB:nearlySingularMatrix')
            disp(['Assilimation_failed_in_Profile_', num2str(iprofile)]);
            % msgid = [];
            continue
        end
        
        % profiles
        if length(wosis_layer_obs) > 1 && isnan(soc_mod(1)) == 0
            optimize_profile_soc = interp1(zsoi(1:soil_decom_num, 1), soc_mod, wosis_layer_depth, 'pchip');
            % scatter(current_obsC(:,1), current_obsC(:,2));
            site_by_site_mod(iprofile_hat, :) = soc_mod;
            indi_map_project_mod(iprofile_hat, 1:length(wosis_layer_obs), control_loc) = optimize_profile_soc;
            indi_map_project_obs(iprofile_hat, 1:length(wosis_layer_obs), control_loc) = wosis_layer_obs;
            indi_map_project_depth(iprofile_hat, 1:length(wosis_layer_obs), control_loc) = wosis_layer_depth;
            
        end
        %% component E
        control_loc = 8;
        % parameters values
        para = nn_predict(iprofile_hat, :);
        para(para_E_loc) = para_global_constant(para_E_loc);
        
        % clear warning info
        lastwarn('');
        soc_mod = ...
            fun_var_decom(para, nbedrock, sand_vector, npp_mean, ...
            input_vector_cwd, input_vector_litter1, input_vector_litter2, input_vector_litter3, ...
            altmax_current_profile, altmax_lastyear_profile, soil_temp_profile, soil_water_profile, xio, xin);
        soc_mod = soc_mod(:, 5);
        
        [~, msgid] = lastwarn;
        if strcmp(msgid, 'MATLAB:singularMatrix') || strcmp(msgid, 'MATLAB:nearlySingularMatrix')
            disp(['Assilimation_failed_in_Profile_', num2str(iprofile)]);
            % msgid = [];
            continue
        end
        
        
        % profiles
        if length(wosis_layer_obs) > 1 && isnan(soc_mod(1)) == 0
            optimize_profile_soc = interp1(zsoi(1:soil_decom_num, 1), soc_mod, wosis_layer_depth, 'pchip');
            % scatter(current_obsC(:,1), current_obsC(:,2));
            site_by_site_mod(iprofile_hat, :) = soc_mod;
            indi_map_project_mod(iprofile_hat, 1:length(wosis_layer_obs), control_loc) = optimize_profile_soc;
            indi_map_project_obs(iprofile_hat, 1:length(wosis_layer_obs), control_loc) = wosis_layer_obs;
            indi_map_project_depth(iprofile_hat, 1:length(wosis_layer_obs), control_loc) = wosis_layer_depth;
            
        end
        
        %% No variantion
        control_loc = 9;
        % parameters values
        para = para_global_constant;
        
        % clear warning info
        lastwarn('');
        soc_mod = ...
            fun_var_decom(para, nbedrock, sand_vector, npp_mean, ...
            input_vector_cwd, input_vector_litter1, input_vector_litter2, input_vector_litter3, ...
            altmax_current_profile, altmax_lastyear_profile, soil_temp_profile, soil_water_profile, xio, xin);
        soc_mod = soc_mod(:, 5);
        
        [~, msgid] = lastwarn;
        if strcmp(msgid, 'MATLAB:singularMatrix') || strcmp(msgid, 'MATLAB:nearlySingularMatrix')
            disp(['Assilimation_failed_in_Profile_', num2str(iprofile)]);
            % msgid = [];
            continue
        end
        
        % profiles
        if length(wosis_layer_obs) > 1 && isnan(soc_mod(1)) == 0
            optimize_profile_soc = interp1(zsoi(1:soil_decom_num, 1), soc_mod, wosis_layer_depth, 'pchip');
            % scatter(current_obsC(:,1), current_obsC(:,2));
            site_by_site_mod(iprofile_hat, :) = soc_mod;
            indi_map_project_mod(iprofile_hat, 1:length(wosis_layer_obs), control_loc) = optimize_profile_soc;
            indi_map_project_obs(iprofile_hat, 1:length(wosis_layer_obs), control_loc) = wosis_layer_obs;
            indi_map_project_depth(iprofile_hat, 1:length(wosis_layer_obs), control_loc) = wosis_layer_depth;
            
        end
    end
    
    
    %% Evaluation Matric (depth)
    compotent_matric = nan(size(indi_map_project_mod, 3), 4);
    
    
    for icomponent = 1:size(indi_map_project_mod, 3)
        
        
        indi_map_project_mod_1d = reshape(indi_map_project_mod(:, :, icomponent), [length(indi_map_project_mod)*150, 1]);
        indi_map_project_obs_1d = reshape(indi_map_project_obs(:, :, icomponent), [length(indi_map_project_obs)*150, 1]);
        
        indi_map_project_obs_1d(indi_map_project_mod_1d < 0) = nan;
        indi_map_project_mod_1d(indi_map_project_mod_1d < 0) = nan;
        
%         valid_layer = find(isnan(indi_map_project_obs_1d) == 0 & indi_map_project_mod_1d > 1);
%         scatter(log10(indi_map_project_obs_1d(valid_layer)), log10(indi_map_project_mod_1d(valid_layer)), 3, 'filled')
%         hold on
%         plot([1.5, 6], [1.5, 6]);
%         xlim([1.5, 6])
%         ylim([1.5, 6])
        
        corr_info = corr(indi_map_project_mod_1d(isnan(indi_map_project_mod_1d) == 0), indi_map_project_obs_1d(isnan(indi_map_project_obs_1d) == 0));
        lm_info = fitlm(indi_map_project_obs_1d(isnan(indi_map_project_obs_1d) == 0), indi_map_project_mod_1d(isnan(indi_map_project_mod_1d) == 0));
        residual_info = var(log((indi_map_project_mod_1d(isnan(indi_map_project_mod_1d) == 0)./indi_map_project_obs_1d(isnan(indi_map_project_obs_1d) == 0))));
        
        % regression with 1:1 line
        ss_tot = sum((indi_map_project_obs_1d(isnan(indi_map_project_mod_1d) == 0) - mean(indi_map_project_obs_1d(isnan(indi_map_project_mod_1d) == 0))).^(2));
        ss_res = sum((indi_map_project_mod_1d(isnan(indi_map_project_mod_1d) == 0) - indi_map_project_obs_1d(isnan(indi_map_project_obs_1d) == 0)).^(2));
        coeff_efficiency = 1 - ss_res/ss_tot;
        
        mse = sum((indi_map_project_mod_1d/1000 - indi_map_project_obs_1d/1000).^2, 'omitnan')/length(indi_map_project_mod_1d);
        
        compotent_matric(icomponent, 1) = coeff_efficiency;
        compotent_matric(icomponent, 2) = corr_info;
        compotent_matric(icomponent, 3) = lm_info.Rsquared.Ordinary;
        compotent_matric(icomponent, 4) = mse;
        
    end
    
    
    %% save simu results
    
    save_pathway = {[data_path, 'output_data/world_simulation_analyses/compontent_control_spatial_variation_', model_name, '_', time_domain, '_', exp_name, '.mat']};
    save_var_namelist = {'compotent_matric'};
    
    save_var = {compotent_matric};
    fun_save_post_da_simu(save_pathway, save_var_namelist, save_var);
    %         save([data_path, 'output_data/world_simulation_analyses/indi_project_', model_name, '_', time_domain, '_', exp_name, '.mat'], 'indi_project');
    %         save([data_path, 'output_data/world_simulation_analyses/compontent_control_spatial_variation_', model_name, '_', time_domain, '_', exp_name, '.mat'], 'compotent_matric');
    
end

disp('Programme Finished');
