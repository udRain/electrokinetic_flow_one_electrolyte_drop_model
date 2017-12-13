function [Ncm_mx,Ncp_mx] = Mx4ChebOPs( r,D1,D2,M,par_dt,par_PE,par_NP,par_CD )
%MX4CHEBOPS Summary of this function goes here
%   Detailed explanation goes here
N=length(r);
Ncm_mx=zeros(N,N,M);
Ncp_mx=zeros(N,N,M);
I=eye(N);
for j=1:M
    grad2=D2+diag(2./r)*D1-j*(j-1)*diag(1./r.^2);
    Ncm_mx(:,:,j)=1/par_dt*I-par_NP*grad2;
    Ncp_mx(:,:,j)=1/par_dt*I-par_NP*grad2;
    Ncm_mx(N,:,j)=D1(N,:)+par_CD/par_PE*I(N,:);
    if mod(j,2)==0
        Ncm_mx(1,:,j)=D1(1,:)+2/r(1)*I(1,:);
    else
        Ncm_mx(1,:,j)=D1(1,:)+3/r(1)*I(1,:);
    end
    Ncm_mx(1,:,1)=I(1,:);
    Ncp_mx(N,:,j)=D1(N,:)-par_CD/par_PE*I(N,:);
    if mod(j,2)==0
        Ncp_mx(1,:,j)=D1(1,:)+2/r(1)*I(1,:); 
    else
        Ncp_mx(1,:,j)=D1(1,:)+3/r(1)*I(1,:);
    end
    Ncp_mx(1,:,1)=I(1,:);
end
end

