%clear all;
%close all;
function a_out = a_matrix(fcwdl2)
global n_soil_layer

nlevdecomp = n_soil_layer;

%use_vertsoilc = 1;
% creat diagnal matrics
a_ma_vr = diag(-ones(4*nlevdecomp, 1));

fcwdl3 = 1 - fcwdl2;

for j = 1:nlevdecomp
    a_ma_vr((3-1)*nlevdecomp+j,(1-1)*nlevdecomp+j) = fcwdl2;
    a_ma_vr((4-1)*nlevdecomp+j,(1-1)*nlevdecomp+j) = fcwdl3;
end
a_out = a_ma_vr;
end
