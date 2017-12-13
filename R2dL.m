function dL_mx = R2dL( R_mx,M )
%REAL_TO_LEGENDRE Summary of this function goes here
%   Detailed explanation goes here
global dP_mx;
global wts;
[n_rpts,~]=size(R_mx);
dL_mx=zeros(n_rpts,M);
for i=2:M
    for j=1:n_rpts
        temp=0;
        for k=1:16
            temp=temp+R_mx(j,k)*dP_mx(i,k)*wts(k);
        end
        dL_mx(j,i)=temp*(2*i-1)/2/(i-1)/i;            
    end
end         
end

