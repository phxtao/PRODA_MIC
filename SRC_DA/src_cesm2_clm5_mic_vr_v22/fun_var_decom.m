function [soc_layer] = fun_var_decom(kp, nbedrock, sand_vector, npp_mean, ...
    input_vector_cwd, input_vector_litter1, input_vector_litter2, input_vector_litter3, ...
    altmax_current_profile, altmax_lastyear_profile, soil_temp_profile, soil_water_profile, xio, xin)
% define global parameters
% global diffus adv cryo q10 fq10 maxpsi minpsi efolding ...
%     tau4cwd tau4l1 tau4l2 tau4l3 tau4s1 tau4s2 tau4s3 fl1s1 fl2s1 fl3s2 fs2s1 fs3s1 ...
%     fs1s2 fs1s3 fs2s3 fcwdl2 fcwdl3 ins beta4cwd beta4l1 beta4l2 beta4l3 dmax p4cwd ...
%     p4ml p4cl p4ll
is_default = 0; 

use_beta = 1; 
month_num = 12;
normalize_q10_to_century_tfunc = false;

global kelvin_to_celsius
kelvin_to_celsius = 273.15;

global use_vertsoilc npool npool_vr n_soil_layer days_per_year secspday
use_vertsoilc = 1 ; % whether or not use vertical maxing part
npool = 7;  % number of pools if no vertical
npool_vr = 140; % number of pools if vertical
n_soil_layer = 20;  % number of soil layers
days_per_year = 365;
secspday = 24*60*60;
% dt = secspday*30;

global max_altdepth_cryoturbation max_depth_cryoturb
max_altdepth_cryoturbation = 2;
max_depth_cryoturb = 3;

global dz dz_node zisoi zsoi
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
zisoi_0 = 0;
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

%-------------------------------------------
% define parameters
%-------------------------------------------
% ------------vertical transport
% define parameter names
% diffusion (bioturbation) 10^(-4) (m2/yr)
bio = kp(1)*(5*10^(-4) - 3*10^(-5)) + 3*10^(-5);
% cryoturbation 5*10^(-4) (m2/yr)
cryo = kp(2)*(16*10^(-4) - 3*10^(-5)) + 3*10^(-5);
% ------------env. scalar
% Q10 (unitless) 1.5
q10 = kp(3)*(3 - 1.2) + 1.2;
% Q10 when forzen (unitless) 1.5
fq10 = q10;
% parameters used in vertical discretization of carbon inputs 10 (metre) 0.5
efolding = 0.5; % kp(4)*(1 - 0) + 0;
% water scaling factor 1
w_scaling = kp(4)*(5 - 0) + 0;
% ------------max decomposition
% turnover time of CWD (yr) 3.3333
tau4cwd = kp(5)*(6 - 1) + 1;
% tau for metabolic litter (yr) 0.0541
tau4l1 = kp(6)*(0.1 - 0) + 0;
% tau for cellulose litter (yr) 0.2041
tau4l2 = kp(7)*(0.3 - 0.1) + 0.1;
% tau for lignin litter (yr)
tau4l3 = tau4l2;
% tau for DOC (yr) 1.1*10-2
tau4s1 = kp(8)*(0.03 - (0.001)) + (0.001); % 10^(kp(9)*((0) - (-3)) + (-3)); % kp(9)*(0.1 - (0.001)) + (0.001); % 
% tau for MIC (yr) 0.57 for death
tau4s2_death = kp(9)*(2 - (0)) + (0); % 10^(kp(10)*(0 - (-2)) + (-2)); % kp(10)*(1 - (0.01)) + (0.01); % 10^(kp(10)*(0 - (-2)) + (-2));
% tau for MIC (yr) 22 for enz production
tau4s2_enz = kp(10)*(30 - (15)) + (15); % kp(11)*(25 - 0.5) + 0.5;
% tau for ENZ (yr) 0.11
tau4s3 = 10^(kp(13)*((0) - (-3)) + (-3)); % kp(11)*(0.15 - 0.0) + 0.0;
% tau for SOC (yr)  default: 4.6*10^-5 -- 1.1*10^-4
tau4s4 = kp(12)*(3*10^(-4) - 0) + 0; % 10^(kp(13)*((-2) - (-5)) + (-5)); % kp(13)*(10^(-3) - 10^(-6)) + 10^(-6); % 10^(kp(13)*((-2) - (-7)) + (-7));
% ------------michaelis-menten
% (michaelis-menten) concentration DOC (yr) for half max assimlation reaction from DOC to MIC, unit: gc/m3, 4*10^2
mm_const_assim = kp(13)*(3000 - 300) + 300; % 10^(kp(14)*(4 - 0) + 0); % kp(14)*(1000 - 10^1) + 10^1; % 10^(kp(14)*(4 - 0) + 0);
% (michaelis-menten) concentration DOC (yr) for half max decomposition reaction from SOC to MIC, unit: gc/m3, 5*10^(4) -- 6*10^5
mm_const_decom = kp(14)*(10^6 - 10^5) + 10^5; % 10^(kp(15)*(9 - 5) + 5); % 10^(kp(15)*(10 - 5) + 5); % kp(15)*(10^6 - 10^4) + 10^4; % 10^(kp(15)*(9 - 4) + 4);
% ------------transfer fraction intra litter
% fraction from cwd to l2 0.75
fcwdl2 = kp(15)*(1 - 0.5) + 0.5;
% ------------transfer fraction litter -> soil
pl1s1 = kp(16)*(0.1 - 0) + 0;
pl2s1 = kp(17)*(0.3 - 0.05) + 0.05;
pl3s4 = kp(18)*(0.95 - 0.6) + 0.6;
% fraction from l1 to s2 0.45
l1_cue = kp(19)*(0.9 - 0.4) + 0.4;
l2_cue = kp(20)*(0.4 - 0) + 0;
l3_cue = l2_cue;

fl1s1 = pl1s1;
fl1s2 = (1-pl1s1)*l1_cue;

fl2s1 = pl2s1;
fl2s2 = (1-pl2s1)*l2_cue;

fl3s2 = (1-pl3s4)*l3_cue;
fl3s4 = pl3s4;
% ------------transfer fraction intra soil
% cue
mic_cue = kp(21)*(0.7 - 0.01) + 0.01;
% fraction from s1 to s2
fs1s2 = mic_cue;
% fraction of cue that leads to death (doc + soc) 0.5
pdeath2soc = kp(22)*(1 - 0) + 0;
% fraction from enz to doc
fs3s1 = 1;
% fraction from soc to doc
fs4s1 = 1;
% ------------input allocation
% beta to describe the shape of vertical profile 0.95
beta = kp(23)*(0.9999 - 0.5) + 0.5;

% maximum and minimum water potential (MPa)
maxpsi= -0.0020;

minpsi= -2;
% parameter for advection (m/yr)
adv = 0;
%% Environmental Scalar (Xi)
xit = nan(n_soil_layer, month_num);

for imonth = 1:month_num
    % temperature related function xit
    % calculate rate constant scalar for soil temperature
    % assuming that the base rate constants are assigned for non-moisture
    % limiting conditions at 25 C.
    for ilayer = 1 : n_soil_layer
        if soil_temp_profile(ilayer, imonth) >= 0 + kelvin_to_celsius
            xit(ilayer, imonth) = q10^((soil_temp_profile(ilayer, imonth) - (kelvin_to_celsius + 25))/10);
        else
            xit(ilayer, imonth) = q10^((273.15 - 298.15)/10)*(fq10^((soil_temp_profile(ilayer, imonth) - (0 + kelvin_to_celsius))/10));
        end
    end
    
    catanf_30 = catanf(30);
    normalization_tref = 15;
    if normalize_q10_to_century_tfunc == true
        % scale all decomposition rates by a constant to compensate for offset between original CENTURY temp func and Q10
        normalization_factor = (catanf(normalization_tref)/catanf_30) / (q10^((normalization_tref-25)/10));
        xit = xit*normalization_factor;
    end
    
end

xiw = soil_water_profile*w_scaling;
xiw(xiw > 1) = 1;


depth_scalar = exp(-zsoi(1:n_soil_layer)/efolding);
monthly_mean_xi_wt = mean(xit, 2).*mean(xiw, 2).*mean(xio, 2).*depth_scalar;
monthly_mean_xin = mean(xin, 2);

%% Triangle Matrix, A Matrix and K Matrix
sand_vector = mean(sand_vector, 2, 'omitnan');
% Allocation matrix
a_ma = a_matrix(fcwdl2);
% decomposition matrix
kk_ma = kk_matrix(monthly_mean_xi_wt, monthly_mean_xin, tau4cwd, tau4l1, tau4l2, tau4l3);

tri_ma_middle = nan(4*20, 4*20, month_num);
for imonth = 1:month_num
    % tri matrix
    monthly_nbedrock = nbedrock(imonth);
    monthly_altmax_current_profile = altmax_current_profile(imonth);
    monthly_altmax_lastyear_profile = altmax_lastyear_profile(imonth);
    tri_ma_middle(:, :, imonth) = tri_matrix(monthly_nbedrock, monthly_altmax_current_profile, monthly_altmax_lastyear_profile, bio, adv, cryo);
end
tri_ma = mean(tri_ma_middle, 3, 'omitnan');


%% Vertical Profile
% in the original beta model in Jackson et al 1996, the unit for the depth
% of the soil is cm (dmax*100)

m_to_cm = 100;
vertical_prof = nan(n_soil_layer, 1);

if altmax_lastyear_profile > 0
    for j = 1:n_soil_layer
        if j == 1
            vertical_prof(j) = (beta^((zisoi_0)*m_to_cm) - beta^(zisoi(j)*m_to_cm))/dz(j);
        else
            vertical_prof(j) = (beta^((zisoi(j - 1))*m_to_cm) - beta^(zisoi(j)*m_to_cm))/dz(j);
        end
    end
else
    vertical_prof(1) = 1/dz(1);
    vertical_prof(2:end) = 0;
end

vertical_input = dz(1:n_soil_layer).*vertical_prof/sum(vertical_prof.*dz(1:n_soil_layer));

%% analytical solution for litter pools
matrix_in=nan(4*n_soil_layer,1);

% calculate annual carbon input, sum monthly input (gc/m3/month) as annual input (gc/m3/year)
input_vector_cwd = sum(input_vector_cwd, 2, 'omitnan');
input_vector_litter1 = sum(input_vector_litter1, 2, 'omitnan');
input_vector_litter2 = sum(input_vector_litter2, 2, 'omitnan');
input_vector_litter3 = sum(input_vector_litter3, 2, 'omitnan');

input_tot_cwd = input_vector_cwd;
input_tot_litter1 = input_vector_litter1;
input_tot_litter2 = input_vector_litter2;
input_tot_litter3 = input_vector_litter3;

matrix_in(1:20,1) = input_tot_cwd*vertical_input./(dz(1:n_soil_layer)*days_per_year); % litter input gc/m3/day
matrix_in(21:40,1) = input_tot_litter1*vertical_input./(dz(1:n_soil_layer)*days_per_year);
matrix_in(41:60,1) = input_tot_litter2*vertical_input./(dz(1:n_soil_layer)*days_per_year);
matrix_in(61:80,1) = input_tot_litter3*vertical_input./(dz(1:n_soil_layer)*days_per_year);

carbon_input = input_tot_cwd + input_tot_litter1 + input_tot_litter2 + input_tot_litter3;

days_in_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]';
site_input_ratio = carbon_input/(sum(npp_mean.*days_in_month*24*3600));

litter_ss = (a_ma*kk_ma-tri_ma)\(-matrix_in); % unit gc m-3
litter_ss(litter_ss < 0) = 0; 

%% analytical solution for soil pools
input_soc = diag(kk_ma(61:80, 61:80)).*litter_ss(61:80)*fl3s4; % unit gc m-3 day-1

input_doc = diag(kk_ma(21:40, 21:40)).*litter_ss(21:40)*fl1s1 + ...
    diag(kk_ma(41:60, 41:60)).*litter_ss(41:60)*fl2s1; % unit gc m-3 day-1

input_mic = diag(kk_ma(21:40, 21:40)).*litter_ss(21:40)*fl1s2 + ...
    diag(kk_ma(41:60, 41:60)).*litter_ss(41:60)*fl2s2 + ...
    diag(kk_ma(61:80, 61:80)).*litter_ss(61:80)*fl3s2; % unit gc m-3 day-1

% microbial organic carbon pool
mic_ss = (input_mic + fs1s2*(input_soc + input_doc))./...
    ((1 - fs1s2).*(1/(tau4s2_death*days_per_year) + 1/(tau4s2_enz*days_per_year)).*monthly_mean_xi_wt);
        
% enzyme carbon pool
enz_ss = (1/(tau4s2_enz*days_per_year))./(1/(tau4s3*days_per_year)).*mic_ss;

% soil organic carbon pool
soc_ss = (input_soc + pdeath2soc * 1/(tau4s2_death*days_per_year) .* monthly_mean_xi_wt .* mic_ss) .* mm_const_decom.*monthly_mean_xi_wt ./ ...
            (1/(tau4s4*days_per_year) .*monthly_mean_xi_wt .* enz_ss - pdeath2soc * 1/(tau4s2_death*days_per_year) .*monthly_mean_xi_wt .* mic_ss - input_soc);

% dissolved organic carbon
doc_ss = ((1/(tau4s2_death*days_per_year) + 1/(tau4s2_enz*days_per_year)).*monthly_mean_xi_wt.*mic_ss.*mm_const_assim.*monthly_mean_xi_wt - input_mic.*mm_const_assim.*monthly_mean_xi_wt) ./ ...
            ((fs1s2 * 1/(tau4s1*days_per_year) .* monthly_mean_xi_wt - (1/(tau4s2_death*days_per_year) + 1/(tau4s2_enz*days_per_year)).*monthly_mean_xi_wt).*mic_ss + input_mic);
  
cpools = [litter_ss(1:20), litter_ss(21:40), litter_ss(41:60), litter_ss(61:80), ...
    doc_ss, mic_ss, enz_ss, soc_ss];

cpools_layer = sum(cpools, 2);
soc_layer = cpools(:, 5:8);
soc_layer = [soc_layer, sum(soc_layer, 2)]; % unit gC/m3

soc_stock = sum(soc_layer.*repmat(dz(1:n_soil_layer), [1, 5]), 1); % unit gC/m2

%% mass balance constranis
% % soc
% mass_balance_index_soc = (1/(tau4s4*days_per_year) .*monthly_mean_xi_wt .* enz_ss)./(mm_const_decom + soc_ss);
% % doc
% mass_balance_index_doc = (1/(tau4s1*days_per_year) .*monthly_mean_xi_wt .* mic_ss)./(mm_const_assim + doc_ss);
% % mic
% mass_balance_index_mic = (1/(tau4s2_death*days_per_year) + 1/(tau4s2_enz*days_per_year)).*monthly_mean_xi_wt;
% % enz
% mass_balance_index_enz = (1/(tau4s3*days_per_year).*monthly_mean_xi_wt);

%% exttra constranis
soil_15cm_loc = 4;


if isempty(find(cpools < 0, 1)) == 0
    
    cpools(cpools > 0) = -cpools(cpools > 0);
    soc_stock(soc_stock > 0) = -soc_stock(soc_stock > 0);
    soc_layer(soc_layer > 0) = -soc_layer(soc_layer > 0);
end

% if isempty(find([mass_balance_index_soc; mass_balance_index_doc; mass_balance_index_mic; mass_balance_index_enz] > 1, 1)) == 0 || ...
%         isempty(find(cpools < 0, 1)) == 0 || ...
%         soc_layer(soil_15cm_loc, 4)/soc_layer(soil_15cm_loc, 2) < 10 || ...
%         soc_layer(soil_15cm_loc, 4)/soc_layer(soil_15cm_loc, 1) < 100 || ...
%         soc_stock(4)/soc_stock(2) < 20 || ... % soc >> mic (total)
%         soc_stock(4)/soc_stock(1) < 10 || ... % soc > doc
%         soc_stock(2)/soc_stock(3) < 1 % mic > enz
%     
%     
%     cpools(cpools > 0) = -cpools(cpools > 0);
%     soc_stock(soc_stock > 0) = -soc_stock(soc_stock > 0);
%     soc_layer(soc_layer > 0) = -soc_layer(soc_layer > 0);
% end

% 
% if isempty(find([mass_balance_index_soc; mass_balance_index_doc; mass_balance_index_mic; mass_balance_index_enz] > 1, 1)) == 0 || ...
%         isempty(find(cpools < 0, 1)) == 0 || ...
%         soc_stock(4)/soc_stock(2) < 10  % soc >> mic (total)    
%     
%     cpools(cpools > 0) = -cpools(cpools > 0);
%     soc_stock(soc_stock > 0) = -soc_stock(soc_stock > 0);
%     soc_layer(soc_layer > 0) = -soc_layer(soc_layer > 0);
% end

% if isempty(find([mass_balance_index_soc; mass_balance_index_doc; mass_balance_index_mic; mass_balance_index_enz] > 1, 1)) == 0
%     cpools(cpools > 0) = -cpools(cpools > 0);
%     soc_stock(soc_stock > 0) = -soc_stock(soc_stock > 0);
%     soc_layer(soc_layer > 0) = -soc_layer(soc_layer > 0);
% end

end

% ----- CENTURY T response function
function catanf_results = catanf(t1)
catanf_results = 11.75 +(29.7 / pi) * atan( pi * 0.031  * ( t1 - 15.4 ));
end







