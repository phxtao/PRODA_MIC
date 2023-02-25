function cost = cost_fun(layer_weight, optimize_profile_soc, wosis_layer_obs)

obs_var = 0.60; % as indicated in WoSIS dataset

cost = nansum((optimize_profile_soc - wosis_layer_obs).^2)./nansum((wosis_layer_obs - mean(wosis_layer_obs)).^2)*5;

% cost = 1*nansum(layer_weight.*((((optimize_profile_soc).^(1) - (wosis_layer_obs).^(1)).^2)...
%     ./(2.*((obs_var*(wosis_layer_obs).^(1)).^2))));

% cost = 1*nansum(layer_weight.*((((optimize_profile_soc./wosis_layer_obs - 1)).^2)));

% cost = 1*nansum(layer_weight.*(((log(optimize_profile_soc + 3) - log(wosis_layer_obs + 3)).^2)...
%     ./(2.*((0.15*log(wosis_layer_obs + 3)).^2))));

end






