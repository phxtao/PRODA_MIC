
function da_result = fun_load_data(da_result_origin, model_name, iprofile, profile_id)
eval(['da_result = da_result_origin.', model_name, '_result_profile_', num2str(iprofile), '_', num2str(profile_id), ';'])
end