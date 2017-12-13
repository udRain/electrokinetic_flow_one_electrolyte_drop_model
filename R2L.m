function L_mx = R2L( R_mx,M )
%REAL_TO_LEGENDRE Summary of this function goes here
%   Detailed explanation goes here
global P_mx;
global wts;
[n_rpts,~]=size(R_mx);
L_mx=zeros(n_rpts,M);
for i=1:M
    for j=1:n_rpts
        temp=0;
        for k=1:16
            temp=temp+R_mx(j,k)*P_mx(i,k)*wts(k);
        end
        L_mx(j,i)=temp*(2*i-1)/2;            
    end
end      
end

