function resubmit_count = fun_resubmit_profiles(start_screening, end_screening, num_profile, model_name)
% num_profile = 111380;
%
% start_screening = 1;
% end_screening = 1440*40;

% model_name = 'clm_cen_vr';

if end_screening > num_profile
    end_screening = num_profile;
end

%% paths
% mac
% data_path = '/Users/phoenix/Google_Drive/Tsinghua_Luo/Projects/DATAHUB/ENSEMBLE/OUTPUT_DATA/';

% server
data_path = ['/GFPS8p/cess11/taof/datahub/ensemble/output_data/mcmc_summary_', model_name, '/'];


para_mean = load([data_path, model_name, '_para_mean.mat'], 'para_mean');
para_mean = para_mean.para_mean;

para_mean = para_mean(start_screening:end_screening, :);

profile_count = (start_screening:end_screening);

resubmit_count = profile_count(isnan(para_mean(:, 1)) == 1);


% stat_r2 = load([data_path, model_name, '_stat_r2.mat'], 'stat_r2');
% stat_r2 = stat_r2.stat_r2;
% 
% stat_r2 = stat_r2(start_screening:end_screening, :);
% 
% profile_count = (start_screening:end_screening);
% 
% resubmit_count = profile_count(stat_r2(:, 1) < 0);

end



