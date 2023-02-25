%=====================================================================
% Calculate the tri_matrix
%=====================================================================
function tri_ma = tri_matrix(nbedrock, altmax, altmax_lastyear, som_diffus, som_adv_flux, cryoturb_diffusion_k)
global dz dz_node zisoi zsoi days_per_year secspday
global max_altdepth_cryoturbation max_depth_cryoturb
global npool n_soil_layer npool_vr

nlevdecomp = n_soil_layer;

epsilon = 1.e-30;

% change the unit from m2/yr to m2/day, in consistant with the unit of input
som_diffus = som_diffus/days_per_year;
som_adv_flux = som_adv_flux/days_per_year;
cryoturb_diffusion_k = cryoturb_diffusion_k/days_per_year;

tri_ma = zeros(4*n_soil_layer, 4*n_soil_layer);

som_adv_coef = zeros(nlevdecomp+1, 1);             	% SOM advective flux (m/day)
som_diffus_coef = zeros(nlevdecomp+1, 1);               % SOM diffusivity due to bio/cryo-turbation (m2/day)
diffus = zeros(nlevdecomp+1, 1);      % diffusivity (m2/day)  (includes spinup correction, if any)
adv_flux = zeros(nlevdecomp+1, 1);     	% advective flux (m/day)  (includes spinup correction, if any)

a_tri_e = zeros(nlevdecomp, 1);        % "a" vector for tridiagonal matrix
b_tri_e = zeros(nlevdecomp, 1);        % "b" vector for tridiagonal matrix
c_tri_e = zeros(nlevdecomp, 1);        % "c" vector for tridiagonal matrix
r_tri_e = zeros(nlevdecomp, 1);        % "r" vector for tridiagonal solution

a_tri_dz = zeros(nlevdecomp, 1);        % "a" vector for tridiagonal matrix with considering the depth
b_tri_dz = zeros(nlevdecomp, 1);        % "b" vector for tridiagonal matrix with considering the depth
c_tri_dz = zeros(nlevdecomp, 1);        % "c" vector for tridiagonal matrix with considering the depth

d_p1_zp1 = zeros(nlevdecomp+1, 1);     % diffusivity/delta_z for next j
% (set to zero for no diffusion)
d_m1_zm1 = zeros(nlevdecomp+1, 1);     % diffusivity/delta_z for previous j
% (set to zero for no diffusion)
f_p1 = zeros(nlevdecomp+1, 1);     % water flux for next j
f_m1 = zeros(nlevdecomp+1, 1);     % water flux for previous j
pe_p1 = zeros(nlevdecomp+1, 1);   % Peclet # for next j
pe_m1 = zeros(nlevdecomp+1, 1);    % Peclet # for previous j

%------ first get diffusivity / advection terms -------%
% use different mixing rates for bioturbation and cryoturbation, with fixed bioturbation and cryoturbation set to a maximum depth
if  (( max(altmax, altmax_lastyear) <= max_altdepth_cryoturbation ) && ...
        ( max(altmax, altmax_lastyear) > 0.) )
    % use mixing profile modified slightly from Koven et al. (2009): constant through active layer, linear decrease from base of active layer to zero at a fixed depth
    for j = 1 : nlevdecomp+1
        if ( j <= nbedrock+1 )
            if ( zisoi(j) < max(altmax, altmax_lastyear) )
                som_diffus_coef(j) = cryoturb_diffusion_k;
                som_adv_coef(j) = 0.;
            else
                som_diffus_coef(j) = max(cryoturb_diffusion_k * ...
                    ( 1. - ( zisoi(j) - max(altmax, altmax_lastyear) ) / ...
                    ( min(max_depth_cryoturb, zisoi(nbedrock+1)) - max(altmax, altmax_lastyear) ) ), 0.);  % go linearly to zero between ALT and max_depth_cryoturb
                som_adv_coef(j) = 0.;
            end
        else
            som_adv_coef(j) = 0.;
            som_diffus_coef(j) = 0.;
        end
    end
elseif (  max(altmax, altmax_lastyear) > 0. )
    % constant advection, constant diffusion
    for j = 1 : nlevdecomp+1
        if ( j <= nbedrock+1 )
            som_adv_coef(j) = som_adv_flux;
            som_diffus_coef(j) = som_diffus;
        else
            som_adv_coef(j) = 0.;
            som_diffus_coef(j) = 0.;
        end
    end
else
    % completely frozen soils--no mixing
    for j = 1 : nlevdecomp+1
        som_adv_coef(j) = 0.;
        som_diffus_coef(j) = 0.;
    end
end

for j = 1 : nlevdecomp+1
    if ( abs(som_adv_coef(j))  < epsilon )
        adv_flux(j) = epsilon;
    else
        adv_flux(j) = som_adv_coef(j);
    end
    if ( abs(som_diffus_coef(j))  < epsilon )
        diffus(j) = epsilon;
    else
        diffus(j) = som_diffus_coef(j);
    end
    
end
% Set Pe (Peclet #) and D/dz throughout column
for j = 1 : nlevdecomp+1
    if ( abs(som_adv_coef(j))  < epsilon )
        adv_flux(j) = epsilon;
    else
        adv_flux(j) = som_adv_coef(j);
    end
    %
    if ( abs(som_diffus_coef(j))  < epsilon )
        diffus(j) = epsilon;
    else
        diffus(j) = som_diffus_coef(j);
    end
    
    
    % Calculate the D and F terms in the Patankar algorithm
    if (j == 1)
        d_m1_zm1(j) = 0.;
        w_p1 = (zsoi(j+1) - zisoi(j)) / dz_node(j+1);
        if ( diffus(j+1) > 0. && diffus(j) > 0.)
            d_p1 = 1. / ((1. - w_p1) / diffus(j) + w_p1 / diffus(j+1)); % Harmonic mean of diffus
        else
            d_p1 = 0.;
        end
        d_p1_zp1(j) = d_p1 / dz_node(j+1);
        f_m1(j) = adv_flux(j);  % Include infiltration here
        f_p1(j) = adv_flux(j+1);
        pe_m1(j) = 0.;
        pe_p1(j) = f_p1(j) / d_p1_zp1(j); % Peclet #
    elseif (j >= nbedrock+1)
        % At the bottom, assume no gradient in d_z (i.e., they're the same)
        w_m1 = (zisoi(j-1) - zsoi(j-1)) / dz_node(j);
        if ( diffus(j) > 0. && diffus(j-1) > 0.)
            d_m1 = 1. / ((1. - w_m1) / diffus(j) + w_m1 / diffus(j-1)); % Harmonic mean of diffus
        else
            d_m1 = 0.;
        end
        d_m1_zm1(j) = d_m1 / dz_node(j);
        d_p1_zp1(j) = d_m1_zm1(j); % Set to be the same
        f_m1(j) = adv_flux(j);
        %f_p1(j) = adv_flux(j+1)
        f_p1(j) = 0.;
        pe_m1(j) = f_m1(j) / d_m1_zm1(j); % Peclet #
        pe_p1(j) = f_p1(j) / d_p1_zp1(j); % Peclet #
    else
        % Use distance from j-1 node to interface with j divided by distance between nodes
        w_m1 = (zisoi(j-1) - zsoi(j-1)) / dz_node(j);
        if ( diffus(j-1) > 0. && diffus(j) > 0.)
            d_m1 = 1. / ((1. - w_m1) / diffus(j) + w_m1 / diffus(j-1)); % Harmonic mean of diffus
        else
            d_m1 = 0.;
        end
        w_p1 = (zsoi(j+1) - zisoi(j)) / dz_node(j+1);
        if ( diffus(j+1) > 0. && diffus(j) > 0.)
            d_p1 = 1. / ((1. - w_p1) / diffus(j) + w_p1 / diffus(j+1)); % Harmonic mean of diffus
        else
            d_p1 = (1. - w_m1) * diffus(j) + w_p1 * diffus(j+1); % Arithmetic mean of diffus
        end
        d_m1_zm1(j) = d_m1 / dz_node(j);
        d_p1_zp1(j) = d_p1 / dz_node(j+1);
        f_m1(j) = adv_flux(j);
        f_p1(j) = adv_flux(j+1);
        pe_m1(j) = f_m1(j) / d_m1_zm1(j); % Peclet #
        pe_p1(j) = f_p1(j) / d_p1_zp1(j); % Peclet #
        
    end
end % j; nlevdecomp

% Calculate the tridiagonal coefficients
for j = 0 : nlevdecomp+1
    
    
    if (j == 0) % top layer (atmosphere)
        %a_tri(j) = 0.
        %b_tri(j) = 1.
        %c_tri(j) = -1.
        %b_tri_e(j) = b_tri(j)
    elseif (j == 1)
        % Set statement functions
        aaa = max (0., (1. - 0.1 * abs(pe_m1(j)))^5);  % A function from Patankar, Table 5.2, pg 95
        a_tri_e(j) = -(d_m1_zm1(j) * aaa + max( f_m1(j), 0.)); % Eqn 5.47 Patankar
        aaa = max (0., (1. - 0.1 * abs(pe_p1(j)))^5);  % A function from Patankar, Table 5.2, pg 95
        c_tri_e(j) = -(d_p1_zp1(j) * aaa + max(-f_p1(j), 0.));
        b_tri_e(j) = -a_tri_e(j) - c_tri_e(j);
    elseif (j < nlevdecomp+1)
        aaa = max (0., (1. - 0.1 * abs(pe_m1(j)))^5); % A function from Patankar, Table 5.2, pg 95
        a_tri_e(j) = -(d_m1_zm1(j) * aaa + max( f_m1(j), 0.)); % Eqn 5.47 Patankar
        aaa = max (0., (1. - 0.1 * abs(pe_p1(j)))^5);  % A function from Patankar, Table 5.2, pg 95
        c_tri_e(j) = -(d_p1_zp1(j) * aaa + max(-f_p1(j), 0.));
        b_tri_e(j) = -a_tri_e(j) - c_tri_e(j);
        
    else % j==nlevdecomp+1; 0 concentration gradient at bottom
        %a_tri(j) = -1.
        %b_tri(j) = 1.
        %c_tri(j) = 0.
    end
end % j; nlevdecomp

for j = 1 : nlevdecomp
    % elements for vertical matrix match unit g/m3,  Tri/dz
    a_tri_dz(j) = a_tri_e(j) / dz(j);
    b_tri_dz(j) = b_tri_e(j) / dz(j);
    c_tri_dz(j) = c_tri_e(j) / dz(j);
    
end



% % no vertical transportation in CWD
% for i = 2 : npool
%     for j = 1 : nlevdecomp
%         if (j == 1)   % upper boundary
%             tri_ma(j+(i-1)*nlevdecomp,j+(i-1)*nlevdecomp) = b_tri_dz(j);
%             tri_ma(j+(i-1)*nlevdecomp,j+1+(i-1)*nlevdecomp) = c_tri_dz(1);
%         end
%         
%         if (j < nlevdecomp + 1)
%             if (j <= nbedrock)
%                 tri_ma(j+(i-1)*nlevdecomp,j-1+(i-1)*nlevdecomp) = a_tri_dz(j);
%                 tri_ma(j+(i-1)*nlevdecomp,j+(i-1)*nlevdecomp) = b_tri_dz(j);
%                 if (j ~= nlevdecomp)
%                     tri_ma(j+(i-1)*nlevdecomp,j+1+(i-1)*nlevdecomp) = c_tri_dz(j);
%                 end
%             end
%         else
%             if (j == nbedrock && j ~= nlevdecomp && j > 1)
%                 tri_ma(j-1+(i-1)*nlevdecomp,j-1+(i-1)*nlevdecomp) = tri_ma(j-1+(i-1)*nlevdecomp,j-1+(i-1)*nlevdecomp) + a_tri_dz;  
%             end
%         end
%         
%     end
% end

% no vertical transportation in CWD
for i = 2 : 4
    for j = 1 : nlevdecomp
        tri_ma(j+(i-1)*nlevdecomp,j+(i-1)*nlevdecomp) = b_tri_dz(j);

        if (j == 1)   % upper boundary
            tri_ma(1+(i-1)*nlevdecomp,1+(i-1)*nlevdecomp) = -c_tri_dz(1);
        end
        if (j == nlevdecomp)  % bottom boundary
            tri_ma(nlevdecomp+(i-1)*nlevdecomp,nlevdecomp+(i-1)*nlevdecomp) = -a_tri_dz(nlevdecomp);
        end

        if (j < nlevdecomp) % avoid tranfer from for example, litr3_20th layer to soil1_1st layer
            tri_ma(j+(i-1)*nlevdecomp,j+1+(i-1)*nlevdecomp) = c_tri_dz(j);
        end

        if (j > 1) % avoid tranfer from for example,soil1_1st layer to litr3_20th layer
            tri_ma(j+(i-1)*nlevdecomp,j-1+(i-1)*nlevdecomp) = a_tri_dz(j);
        end
    end
end
end
