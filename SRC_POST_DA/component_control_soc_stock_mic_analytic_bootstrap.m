% clear;
% clc;
% start_id = 1;
% end_id = 1;

%% cesm2 Settings
% parallel setting
delete(gcp('nocreate'));
parpool(10);

cesm2_case_name = 'sasu_f05_g16_checked_step4';
start_year = 661;
end_year = 680;

model_name = 'cesm2_clm5_mic_vr_v22';

nn_exp_name = 'exp_pc_cesm2_24';

time_domain = 'whole_time'; % 'whole_time', 'before_1985', 'after_1985'

% for ibootstrap = start_id:10:end_id % cross_valid_id
parfor ibootstrap = start_id:end_id % cross_valid_id
    disp(['Processing ibootstrap: ', num2str(ibootstrap)]);
    
    nn_exp_name = ['exp_pc_cesm2_24_bootstrap_', num2str(ibootstrap)];
    % nn_exp_name = ['exp_pc_cesm2_24_cross_valid_0_', num2str(ibootstrap)];
    
    % paths
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
    
    
    warning('off');
    format long e;
    
    %% Load NN Prediction Results
    GlobalGrid = load([data_path, 'input_data/data4nn/world_grid_envinfo_present.mat']);
    GlobalGrid = GlobalGrid.EnvInfo;
    
    grid_para_result = csvread([data_path, 'output_data/neural_networking/grid_para_result_', model_name, '_', time_domain, '_', nn_exp_name, '.csv']);
    valid_grid_loc = csvread([data_path, 'output_data/neural_networking/valid_grid_loc_', model_name, '_', time_domain, '_', nn_exp_name, '.csv']);
    
    %% load MCMC MLE results
    para_global_mean = load([data_path, 'output_data/mcmc_summary_', model_name, '/', model_name , '_mcmc_results.mat']);
    para_global_mean = para_global_mean.mcmc_results.para_mean;
    para_global_mean = real(para_global_mean);
    
    para_global_constant = mean(para_global_mean, 1, 'omitnan')';
    % grid info
    GlobalGrid = GlobalGrid(valid_grid_loc, :);
    GridInfo = GlobalGrid(:, [1:2, 12]);
    ValidGrid = length(GridInfo(:, 1));
    
    Var_Decom_Mean_30cm = nan(ValidGrid, 9);
    Var_Decom_Mean_100cm = nan(ValidGrid, 9);
    Var_Decom_Mean_200cm = nan(ValidGrid, 9);
    Var_Decom_Mean_Total = nan(ValidGrid, 9);
    
    online_soc_stock = nan(ValidGrid, 3);
    cesm2_simu_soc_stock =mean(cesm2_simu_soc_stock, 3, 'omitnan');
    
    % cesm2 resolution
    cesm2_resolution_lat = 180/384;
    cesm2_resolution_lon = 360/576;
    LonGrid_Range = (-180 + cesm2_resolution_lon/2 : cesm2_resolution_lon : 180 - cesm2_resolution_lon/2)';
    LatGrid_Range = (90 - cesm2_resolution_lat/2 : -cesm2_resolution_lat : -90 + cesm2_resolution_lat/2)';
    
    for iGrid = 1:length(GridInfo)
        
        warning off
        
        % disp(['Processing Grid: ', num2str(iGrid)]);
        
        Grid_Lon = GridInfo(iGrid,1);
        Grid_Lat = GridInfo(iGrid,2);
        
        lat_loc = find(abs(Grid_Lat - LatGrid_Range) == min(abs(Grid_Lat - LatGrid_Range)));
        lon_loc = find(abs(Grid_Lon - LonGrid_Range) == min(abs(Grid_Lon - LonGrid_Range)));
        
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
            disp(['No input info ', ' Grid ', num2str(iGrid)]);
            continue
        end
        
        % npp from CESM2 simulatoin
        npp_mean = reshape(cesm2_simu_npp(lat_loc, lon_loc, :), [month_num, 1]);
        % NPP info of studied profile from soc_modIS NPP mean (year 2001-2016)
        % npp_mean = env_info(iprofile, 15);
        
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
        
        % Parameter names and their initial values in MCMC
        para_name = {'bio', 'cryo', 'q10', 'w_scaling', ...
            'tau4cwd', 'tau4l1', 'tau4l2', 'tau4s1', 'tau4s2_death', 'tau4s2_enz', 'tau4s3', 'tau4s4', ...
            'mm_const_assim', 'mm_const_decom', ...
            'fcwdl2', ...
            'pl1s1', 'pl2s1', 'pl3s4', 'l1_cue', 'l2l3_cue', ...
            'mic_cue', 'pdeath2soc', ...
            'beta'};
        
        para_A_loc = [8, 13, 19:21]; 
        para_K_loc = [5:7, 9:12, 14]; 
        para_Xi_loc = [3, 4];
        para_I_loc = 23;
        para_V_loc = [1, 2];
        para_E_loc = [15:18, 22]; 
        
        %% fully varying variable
        Para = grid_para_result(iGrid, :);
        
        soc_mod = fun_var_decom(Para, nbedrock, sand_vector, npp_mean, ...
            input_vector_cwd, input_vector_litter1, input_vector_litter2, input_vector_litter3, ...
            altmax_current_profile, altmax_lastyear_profile, soil_temp_profile, soil_water_profile, xio, xin);
        soc_layer = soc_mod(:, 5);
        
        var_loc = 1;
        
        Var_Decom_Mean_30cm(iGrid, var_loc) = sum(soc_layer(1:4).*dz(1:4)) + soc_layer(5)*dz(5)*(0.3 - zisoi(4))/(zisoi(5) - zisoi(4));
        Var_Decom_Mean_100cm(iGrid, var_loc) = sum(soc_layer(1:8).*dz(1:8)) + soc_layer(9)*dz(9)*(1 - zisoi(8))/(zisoi(9) - zisoi(8));
        Var_Decom_Mean_200cm(iGrid, var_loc) = sum(soc_layer(1:11).*dz(1:11)) + soc_layer(12)*dz(12)*(2 - zisoi(11))/(zisoi(12) - zisoi(11));
        Var_Decom_Mean_Total(iGrid, var_loc) = sum(soc_layer(1:20).*dz(1:20));
        
        %% component A: parameter 8 - 16
        Para = grid_para_result(iGrid, :);
        Para(para_A_loc) = para_global_constant(para_A_loc);
        
        soc_mod = fun_var_decom(Para, nbedrock, sand_vector, npp_mean, ...
            input_vector_cwd, input_vector_litter1, input_vector_litter2, input_vector_litter3, ...
            altmax_current_profile, altmax_lastyear_profile, soil_temp_profile, soil_water_profile, xio, xin);
        soc_layer = soc_mod(:, 5);
        
        var_loc = 2;
        
        Var_Decom_Mean_30cm(iGrid, var_loc) = sum(soc_layer(1:4).*dz(1:4)) + soc_layer(5)*dz(5)*(0.3 - zisoi(4))/(zisoi(5) - zisoi(4));
        Var_Decom_Mean_100cm(iGrid, var_loc) = sum(soc_layer(1:8).*dz(1:8)) + soc_layer(9)*dz(9)*(1 - zisoi(8))/(zisoi(9) - zisoi(8));
        Var_Decom_Mean_200cm(iGrid, var_loc) = sum(soc_layer(1:11).*dz(1:11)) + soc_layer(12)*dz(12)*(2 - zisoi(11))/(zisoi(12) - zisoi(11));
        Var_Decom_Mean_Total(iGrid, var_loc) = sum(soc_layer(1:20).*dz(1:20));

        %% component I: parameter 17
        Para = grid_para_result(iGrid, :);
        Para(para_I_loc) = para_global_constant(para_I_loc);
        
        soc_mod = fun_var_decom(Para, nbedrock, sand_vector, npp_mean, ...
            input_vector_cwd, input_vector_litter1, input_vector_litter2, input_vector_litter3, ...
            altmax_current_profile, altmax_lastyear_profile, soil_temp_profile, soil_water_profile, xio, xin);
        soc_layer = soc_mod(:, 5);
        
        var_loc = 3;
        
        Var_Decom_Mean_30cm(iGrid, var_loc) = sum(soc_layer(1:4).*dz(1:4)) + soc_layer(5)*dz(5)*(0.3 - zisoi(4))/(zisoi(5) - zisoi(4));
        Var_Decom_Mean_100cm(iGrid, var_loc) = sum(soc_layer(1:8).*dz(1:8)) + soc_layer(9)*dz(9)*(1 - zisoi(8))/(zisoi(9) - zisoi(8));
        Var_Decom_Mean_200cm(iGrid, var_loc) = sum(soc_layer(1:11).*dz(1:11)) + soc_layer(12)*dz(12)*(2 - zisoi(11))/(zisoi(12) - zisoi(11));
        Var_Decom_Mean_Total(iGrid, var_loc) = sum(soc_layer(1:20).*dz(1:20));
        
        
        %% component K: parameter 5 - 7
        Para = grid_para_result(iGrid, :);
        Para(para_K_loc) = para_global_constant(para_K_loc);
        
        soc_mod = fun_var_decom(Para, nbedrock, sand_vector, npp_mean, ...
            input_vector_cwd, input_vector_litter1, input_vector_litter2, input_vector_litter3, ...
            altmax_current_profile, altmax_lastyear_profile, soil_temp_profile, soil_water_profile, xio, xin);
        soc_layer = soc_mod(:, 5);
        
        
        var_loc = 4;
        
        Var_Decom_Mean_30cm(iGrid, var_loc) = sum(soc_layer(1:4).*dz(1:4)) + soc_layer(5)*dz(5)*(0.3 - zisoi(4))/(zisoi(5) - zisoi(4));
        Var_Decom_Mean_100cm(iGrid, var_loc) = sum(soc_layer(1:8).*dz(1:8)) + soc_layer(9)*dz(9)*(1 - zisoi(8))/(zisoi(9) - zisoi(8));
        Var_Decom_Mean_200cm(iGrid, var_loc) = sum(soc_layer(1:11).*dz(1:11)) + soc_layer(12)*dz(12)*(2 - zisoi(11))/(zisoi(12) - zisoi(11));
        Var_Decom_Mean_Total(iGrid, var_loc) = sum(soc_layer(1:20).*dz(1:20));
        
        
        %% component V: parameter 1 - 2
        Para = grid_para_result(iGrid, :);
        Para(para_V_loc) = para_global_constant(para_V_loc);
        
        soc_mod = fun_var_decom(Para, nbedrock, sand_vector, npp_mean, ...
            input_vector_cwd, input_vector_litter1, input_vector_litter2, input_vector_litter3, ...
            altmax_current_profile, altmax_lastyear_profile, soil_temp_profile, soil_water_profile, xio, xin);
        soc_layer = soc_mod(:, 5);
        
        
        var_loc = 5;
        
        Var_Decom_Mean_30cm(iGrid, var_loc) = sum(soc_layer(1:4).*dz(1:4)) + soc_layer(5)*dz(5)*(0.3 - zisoi(4))/(zisoi(5) - zisoi(4));
        Var_Decom_Mean_100cm(iGrid, var_loc) = sum(soc_layer(1:8).*dz(1:8)) + soc_layer(9)*dz(9)*(1 - zisoi(8))/(zisoi(9) - zisoi(8));
        Var_Decom_Mean_200cm(iGrid, var_loc) = sum(soc_layer(1:11).*dz(1:11)) + soc_layer(12)*dz(12)*(2 - zisoi(11))/(zisoi(12) - zisoi(11));
        Var_Decom_Mean_Total(iGrid, var_loc) = sum(soc_layer(1:20).*dz(1:20));
        
        %% component Xi: parameter 3, 4, and 18
        Para = grid_para_result(iGrid, :);
        Para(para_Xi_loc) = para_global_constant(para_Xi_loc);
        
        soc_mod = fun_var_decom(Para, nbedrock, sand_vector, npp_mean, ...
            input_vector_cwd, input_vector_litter1, input_vector_litter2, input_vector_litter3, ...
            altmax_current_profile, altmax_lastyear_profile, soil_temp_profile, soil_water_profile, xio, xin);
        soc_layer = soc_mod(:, 5);
        
        
        var_loc = 6;
        
        Var_Decom_Mean_30cm(iGrid, var_loc) = sum(soc_layer(1:4).*dz(1:4)) + soc_layer(5)*dz(5)*(0.3 - zisoi(4))/(zisoi(5) - zisoi(4));
        Var_Decom_Mean_100cm(iGrid, var_loc) = sum(soc_layer(1:8).*dz(1:8)) + soc_layer(9)*dz(9)*(1 - zisoi(8))/(zisoi(9) - zisoi(8));
        Var_Decom_Mean_200cm(iGrid, var_loc) = sum(soc_layer(1:11).*dz(1:11)) + soc_layer(12)*dz(12)*(2 - zisoi(11))/(zisoi(12) - zisoi(11));
        Var_Decom_Mean_Total(iGrid, var_loc) = sum(soc_layer(1:20).*dz(1:20));
        
        %% component NPP
        Para = grid_para_result(iGrid, :);
        
        while true
            npp_scale_factor = 1 + randi([-1, 1])*2*(GlobalGrid(iGrid, 73)/GlobalGrid(iGrid, 72));
            if npp_scale_factor ~= 1
                break
            end
        end
        
        soc_mod = fun_var_decom(Para, nbedrock, sand_vector, npp_mean*npp_scale_factor, ...
            input_vector_cwd*npp_scale_factor, input_vector_litter1*npp_scale_factor, input_vector_litter2*npp_scale_factor, input_vector_litter3*npp_scale_factor, ...
            altmax_current_profile, altmax_lastyear_profile, soil_temp_profile, soil_water_profile, xio, xin);
        soc_layer = soc_mod(:, 5);
        
        
        var_loc = 7;
        
        Var_Decom_Mean_30cm(iGrid, var_loc) = sum(soc_layer(1:4).*dz(1:4)) + soc_layer(5)*dz(5)*(0.3 - zisoi(4))/(zisoi(5) - zisoi(4));
        Var_Decom_Mean_100cm(iGrid, var_loc) = sum(soc_layer(1:8).*dz(1:8)) + soc_layer(9)*dz(9)*(1 - zisoi(8))/(zisoi(9) - zisoi(8));
        Var_Decom_Mean_200cm(iGrid, var_loc) = sum(soc_layer(1:11).*dz(1:11)) + soc_layer(12)*dz(12)*(2 - zisoi(11))/(zisoi(12) - zisoi(11));
        Var_Decom_Mean_Total(iGrid, var_loc) = sum(soc_layer(1:20).*dz(1:20));
        
        %% component E
        Para = grid_para_result(iGrid, :);
        Para(para_E_loc) = para_global_constant(para_E_loc);
        
        soc_mod = fun_var_decom(Para, nbedrock, sand_vector, npp_mean, ...
            input_vector_cwd, input_vector_litter1, input_vector_litter2, input_vector_litter3, ...
            altmax_current_profile, altmax_lastyear_profile, soil_temp_profile, soil_water_profile, xio, xin);
        soc_layer = soc_mod(:, 5);
        
        var_loc = 8;
        
        Var_Decom_Mean_30cm(iGrid, var_loc) = sum(soc_layer(1:4).*dz(1:4)) + soc_layer(5)*dz(5)*(0.3 - zisoi(4))/(zisoi(5) - zisoi(4));
        Var_Decom_Mean_100cm(iGrid, var_loc) = sum(soc_layer(1:8).*dz(1:8)) + soc_layer(9)*dz(9)*(1 - zisoi(8))/(zisoi(9) - zisoi(8));
        Var_Decom_Mean_200cm(iGrid, var_loc) = sum(soc_layer(1:11).*dz(1:11)) + soc_layer(12)*dz(12)*(2 - zisoi(11))/(zisoi(12) - zisoi(11));
        Var_Decom_Mean_Total(iGrid, var_loc) = sum(soc_layer(1:20).*dz(1:20));
        
        %% No variantion
        Para = para_global_constant;
        
        soc_mod = fun_var_decom(Para, nbedrock, sand_vector, npp_mean, ...
            input_vector_cwd, input_vector_litter1, input_vector_litter2, input_vector_litter3, ...
            altmax_current_profile, altmax_lastyear_profile, soil_temp_profile, soil_water_profile, xio, xin);
        soc_layer = soc_mod(:, 5);
        
        
        var_loc = 9;
        
        Var_Decom_Mean_30cm(iGrid, var_loc) = sum(soc_layer(1:4).*dz(1:4)) + soc_layer(5)*dz(5)*(0.3 - zisoi(4))/(zisoi(5) - zisoi(4));
        Var_Decom_Mean_100cm(iGrid, var_loc) = sum(soc_layer(1:8).*dz(1:8)) + soc_layer(9)*dz(9)*(1 - zisoi(8))/(zisoi(9) - zisoi(8));
        Var_Decom_Mean_200cm(iGrid, var_loc) = sum(soc_layer(1:11).*dz(1:11)) + soc_layer(12)*dz(12)*(2 - zisoi(11))/(zisoi(12) - zisoi(11));
        Var_Decom_Mean_Total(iGrid, var_loc) = sum(soc_layer(1:20).*dz(1:20));
        
        
    end
    
    component_num = 9;
    depth_num = 4;
    
    % correction by coarse fragment fraction
    raw_soc_0_30cm = Var_Decom_Mean_30cm;
    raw_soc_30_100cm = Var_Decom_Mean_100cm - Var_Decom_Mean_30cm;
    raw_soc_100_200cm = Var_Decom_Mean_200cm - Var_Decom_Mean_100cm;
    raw_soc_200cm_max = Var_Decom_Mean_Total - Var_Decom_Mean_200cm;
    
    for icomponent = 1:component_num
        Var_Decom_Mean_30cm(:, icomponent) = raw_soc_0_30cm(:, icomponent).*(1 - GlobalGrid(:, 45)/100); % coarse fragment fraction for 0 - 30cm
        Var_Decom_Mean_100cm(:, icomponent) = raw_soc_0_30cm(:, icomponent).*(1 - GlobalGrid(:, 45)/100) + ...
            raw_soc_30_100cm(:, icomponent).*(1 - GlobalGrid(:, 46)/100); % coarse fragment fraction for 30 - 100cm
        Var_Decom_Mean_200cm(:, icomponent) = raw_soc_0_30cm(:, icomponent).*(1 - GlobalGrid(:, 45)/100) + ...
            raw_soc_30_100cm(:, icomponent).*(1 - GlobalGrid(:, 46)/100) + ...
            raw_soc_100_200cm(:, icomponent).*(1 - GlobalGrid(:, 47)/100); % coarse fragment fraction for 100 - 200cm
        Var_Decom_Mean_Total(:, icomponent) = raw_soc_0_30cm(:, icomponent).*(1 - GlobalGrid(:, 45)/100) + ...
            raw_soc_30_100cm(:, icomponent).*(1 - GlobalGrid(:, 46)/100) + ...
            raw_soc_100_200cm(:, icomponent).*(1 - GlobalGrid(:, 47)/100) + ...
            raw_soc_200cm_max(:, icomponent).*(1 - GlobalGrid(:, 47)/100); % coarse fragment fraction for 100 - 200cm
    end
    
    global_stock = nan(depth_num, component_num);
    
    Radius = 6371008.8;
    Resolution_Lon = 0.5;
    Resolution_Lat = 0.5;
    
    Length_Top = (2*pi*Radius*cos(abs(GridInfo(:, 2)+Resolution_Lat/2)/180*pi)/360)*Resolution_Lon;
    Length_Down = (2*pi*Radius*cos(abs(GridInfo(:, 2)-Resolution_Lat/2)/180*pi)/360)*Resolution_Lon;
    Height = (pi*Radius/180)*Resolution_Lat;
    grid_area = (Length_Top + Length_Down)*Height/2;
    
    valid_loc = find(Var_Decom_Mean_Total(:, 1) > 0 & Var_Decom_Mean_Total(:, 1) < 1000000);
    
    for icomponent = 1:component_num
        if icomponent == 1
            grid_stock_30cm = Var_Decom_Mean_30cm(valid_loc, icomponent).*grid_area(valid_loc); % unit gC
            grid_stock_100cm = Var_Decom_Mean_100cm(valid_loc, icomponent).*grid_area(valid_loc);
            grid_stock_200cm = Var_Decom_Mean_200cm(valid_loc, icomponent).*grid_area(valid_loc);
            grid_stock_total = Var_Decom_Mean_Total(valid_loc, icomponent).*grid_area(valid_loc);
            
        else
            grid_stock_30cm = abs(Var_Decom_Mean_30cm(valid_loc, icomponent) - Var_Decom_Mean_30cm(valid_loc, 1)).*grid_area(valid_loc); % unit gC
            grid_stock_100cm = abs(Var_Decom_Mean_100cm(valid_loc, icomponent) - Var_Decom_Mean_100cm(valid_loc, 1)).*grid_area(valid_loc);
            grid_stock_200cm = abs(Var_Decom_Mean_200cm(valid_loc, icomponent) - Var_Decom_Mean_200cm(valid_loc, 1)).*grid_area(valid_loc);
            grid_stock_total = abs(Var_Decom_Mean_Total(valid_loc, icomponent) - Var_Decom_Mean_Total(valid_loc, 1)).*grid_area(valid_loc);
        end
        global_stock(1, icomponent) = sum(grid_stock_30cm, 'omitnan')/(10^15); % unit Pg C
        global_stock(2, icomponent) = sum(grid_stock_100cm, 'omitnan')/(10^15); % unit Pg C
        global_stock(3, icomponent) = sum(grid_stock_200cm, 'omitnan')/(10^15); % unit Pg C
        global_stock(4, icomponent) = sum(grid_stock_total, 'omitnan')/(10^15); % unit Pg C
    end
    
    
    save_pathway = {[data_path, 'output_data/world_simulation_analyses/compontent_control_soc_total_stock_',  model_name, '_', nn_exp_name, '_control_test.mat']};
    save_var_namelist = {'global_stock'};
    
    save_var = {global_stock};
    fun_save_post_da_simu(save_pathway, save_var_namelist, save_var);
    
end

disp('Global Projection Finished');
