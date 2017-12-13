function [u_R,v_R,u_L,v_L]=Calc_u_v( r,psi_L,dpsi_L )
%CALC_U_V Summary of this function goes here
%   Detailed explanation goes here
global P_mx;
global dP_mx;

[n_rpts,n_modes]=size(psi_L);

u_L=zeros(n_rpts,n_modes);
v_L=zeros(n_rpts,n_modes);

%for j=1:n_modes
for j=1:n_modes
    u_L(:,j)=-j*(j-1)*psi_L(:,j)./(r.^2);
end
for i=1:n_rpts
    v_L(i,:)=-dpsi_L(i,:)/r(i);
end

u_R=u_L*P_mx(1:n_modes,:);
v_R=v_L*dP_mx(1:n_modes,:);

end

