function [ r_comp,th_comp ] = Calc_grad_R( r,x_L,D1 )
%CALC_GRADPHI Summary of this function goes here
%   Detailed explanation goes here
global P_mx;
global dP_mx;
dx_dr=D1*x_L;
[~,n_modes]=size(x_L);
r_comp=dx_dr*P_mx(1:n_modes,:);
th_comp=diag(1./r)*x_L*dP_mx(1:n_modes,:);
end

