function     fun_save_da(model_name, iprofile, profile_id, wosis_layer_obs, soc_mod_opt, soc_std_opt, para_summary_mean, para_summary_std, para_es, GR, parameters_keep2, stat_parallel, upgrade2)

eval([model_name, '_result_profile_', num2str(iprofile), '_', num2str(profile_id), '.wosis_layer_obs = wosis_layer_obs;']);

eval([model_name, '_result_profile_', num2str(iprofile), '_', num2str(profile_id), '.para_summary_mean = para_summary_mean;']);
eval([model_name, '_result_profile_', num2str(iprofile), '_', num2str(profile_id), '.para_summary_std = para_summary_std;']);

eval([model_name, '_result_profile_', num2str(iprofile), '_', num2str(profile_id), '.soc_mod_opt = soc_mod_opt;']);
eval([model_name, '_result_profile_', num2str(iprofile), '_', num2str(profile_id), '.soc_std_opt = soc_std_opt;']);

eval([model_name, '_result_profile_', num2str(iprofile), '_', num2str(profile_id), '.para_es = para_es;']);

eval([model_name, '_result_profile_', num2str(iprofile), '_', num2str(profile_id), '.GR = GR;']);

% eval([model_name, '_result_profile_', num2str(iprofile), '_', num2str(profile_id), '.parameters_keep2 = parameters_keep2;']);
eval([model_name, '_result_profile_', num2str(iprofile), '_', num2str(profile_id), '.stat_parallel = stat_parallel;']);

eval([model_name, '_result_profile_', num2str(iprofile), '_', num2str(profile_id), '.upgrade2 = upgrade2;']);

save(['/GFPS8p/cess11/taof/datahub/ensemble/output_data/', model_name, '/',...
    model_name, '_result_profile_', num2str(iprofile), '_', num2str(profile_id), '.mat'],...
    [model_name, '_result_profile_', num2str(iprofile), '_', num2str(profile_id)]);

% save(['/Users/phoenix/Google_Drive/Tsinghua_Luo/Projects/DATAHUB/ENSEMBLE/OUTPUT_FILE/', model_name, '/',...
%     model_name, '_result_profile_', num2str(iprofile), '_', num2str(profile_id), '.mat'],...
%     [model_name, '_result_profile_', num2str(iprofile), '_', num2str(profile_id)]);

end
