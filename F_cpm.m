function [Fres,phi] = F_cpm( r,cm1,cp1,lbc,par_PE,D1,D2 )
%F_CPM Summary of this function goes here
%   Detailed explanation goes here
N=length(r);
I=eye(N);
IN=D1;
IN(1,:)=I(1,:);
phipp=(cm1-cp1)/2/par_PE^2;

phi1=IN\(r.^2.*[0;phipp(2:N)]);
phip=phi1./r.^2-(phi1(N)-lbc)./r.^2;
phi2=IN\(1./r.^2.*[0;phi1(2:N)]);
phi=phi2+0*(phi1(N)-lbc)./r;
%Fcm=cm1-exp(phi);
%Fcp=cp1-exp(-phi);

Fcm=((D2+diag(2./r)*D1)-(diag(phipp)+diag(phip)*D1))*cm1;
Fcp=((D2+diag(2./r)*D1)+(diag(phipp)+diag(phip)*D1))*cp1;
%Fcm(1)=cm2(1)-cp2(1);
%Fcp(1)=D1(1,:)*cp1+1/r(1)*(cp1(1)-1);
%Fcm(1)=D1(1,:)*cm1+1/r(1)*(cm1(1)-1);
Fcm(1)=cm1(1)-1;
Fcp(1)=cp1(1)-1;
%Fcm(N)=D1(N,:)*(cm1+cp1)+lbc*(cp1(N)-cm1(N));
%Fcp(N)=D1(N,:)*(cp1-cm1)+lbc*(cp1(N)+cm1(N));
Fcm(N)=D1(N,:)*cm1-lbc*cm1(N);
Fcp(N)=D1(N,:)*cp1+lbc*cp1(N);

Fres=[Fcm;Fcp];
end

