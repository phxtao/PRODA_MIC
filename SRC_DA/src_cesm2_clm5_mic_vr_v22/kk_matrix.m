%clear all;
%close all;
function kk_out = kk_matrix(monthly_mean_xi_wt, monthly_mean_xin, tau4cwd, tau4l1, tau4l2, tau4l3)
global n_soil_layer days_per_year

kl1 = 1/(days_per_year * tau4l1);
kl2 = 1/(days_per_year  * tau4l2);
kl3 = 1/(days_per_year  * tau4l3);
kcwd = 1/(days_per_year  * tau4cwd);


kk_ma_vr = zeros(4*n_soil_layer);   % kk_matrix, decay matrix * scalar matrix
for j = 1:n_soil_layer
    % CWD exists only on the surface of land
    kk_ma_vr(j,j) = kcwd * monthly_mean_xi_wt(j);
    % other litters and SOC and be decomposed at each layer of the soil
    % profile
    kk_ma_vr(1*n_soil_layer+j, 1*n_soil_layer+j) = kl1 * monthly_mean_xi_wt(j) * monthly_mean_xin(j);
    kk_ma_vr(2*n_soil_layer+j, 2*n_soil_layer+j) = kl2 * monthly_mean_xi_wt(j) * monthly_mean_xin(j);
    kk_ma_vr(3*n_soil_layer+j, 3*n_soil_layer+j) = kl3 * monthly_mean_xi_wt(j) * monthly_mean_xin(j);
end
kk_out = kk_ma_vr;
end
