function setup_modes
%MODES_SETUP Summary of this function goes here
%   Detailed explanation goes here
global P_mx;        %matrix of Legendre polynomial P_j
global dP_mx;       %matrix of 1st derivative of Legendre polynomial dP_j/dtheta
global absc;        %Abscissae of Gauss Legendre quadrature
global wts;         %Weights of Gauss Legendre quadrature
[absc,wts]=Calc_GLQ(16,-1,1,1);
P_mx=zeros(15,16);
dP_mx=zeros(15,16);
for i=1:16
    [P_mx(:,i),dP_mx(:,i)]...
        =Calc_LP(14,absc(i));
end
end

