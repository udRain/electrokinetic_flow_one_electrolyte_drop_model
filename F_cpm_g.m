function [Fres,Lcp_R,Lcm_R,Lcp_L,Lcm_L,phi] = F_cpm_g( r,n_modes,cm1,cp1,cm2,cp2,par_PE,par_CD,par_EP,par_NP,D1,D2 )
%F_CPM_G Summary of this function goes here
%   Detailed explanation goes here

global P_mx;        %matrix of Legendre polynomial P_j
global dP_mx;       %matrix of 1st derivative of Legendre polynomial dP_j/dtheta
global absc;        %Abscissae of Gauss Legendre quadrature
global wts;         %Weights of Gauss Legendre quadrature

N=length(r);
I=eye(N);
IN=D1;
IN(1,:)=I(1,:);

rhspe=1/2*(cm1-cp1);
%{
for k=2:n_modes
    rhspe(:,k)=rhspe(:,k)-rhspe(1,k)*r(1)^k./r.^k;
end
%}
del2phi_L=rhspe/par_PE^2;
del2phi_R=del2phi_L*P_mx(1:n_modes,:);
rhspeb=zeros(1,n_modes);
rhspeb(1)=-par_CD/par_PE;
rhspeb(2)=par_EP;
phi=zeros(N,n_modes);
phip=zeros(N,n_modes);
sl=zeros(N,n_modes);

for j=1:n_modes
    inp1=IN\(r.^(j+1).*[0;del2phi_L(2:N,j)]);
    inp2p=r.^(-2*j).*inp1;
    inp2=IN\([0;inp2p(2:N)]);
    c1_b=(rhspeb(j)-inp1(N)-(j-1)*inp2(N))/j;
    phi(:,j)=r.^(j-1).*inp2-c1_b./r.^j;
    phip(:,j)=(j-1)*r.^(j-2).*inp2+r.^(j-1).*inp2p+j*c1_b./r.^(j+1);
end
phi(:,2)=phi(:,2)-par_EP*r;
phip(:,2)=phip(:,2)-par_EP;
gradphi_r_R=phip*P_mx(1:n_modes,:);
gradphi_th_R=diag(1./r)*phi*dP_mx(1:n_modes,:);
%{
force_r_Lc=zeros(N,n_modes);
force_th_Lc=zeros(N,n_modes);
%force_r_Lc(:,1)=-phip(:,1).*del2phi_L(:,1)-1/3*phip(:,2).*del2phi_L(:,2)-...
%    1/5*phip(:,3).*del2phi_L(:,3);

force_r_Lc(:,2)=-(phip(:,1).*del2phi_L(:,2)+phip(:,2).*del2phi_L(:,1));
force_th_Lc(:,2)=-phi(:,2).*del2phi_L(:,1);
dr_force_th_Lc=zeros(N,n_modes);
dr_force_th_Lc(:,2)=-phip(:,2).*del2phi_L(:,1)-...
    phi(:,2).*drdel2phi_L(:,1).*(sign(r(floor(N/8))-r)+1)/2;
%}
%{
force_r_Lc(:,2)=-(phip(:,1).*del2phi_L(:,2)+phip(:,2).*del2phi_L(:,1))-...
    2/5*(phip(:,3).*del2phi_L(:,2)+phip(:,2).*del2phi_L(:,3));
force_r_Lc(:,3)=-(phip(:,1).*del2phi_L(:,3)+phip(:,3).*del2phi_L(:,1))-...
    2/3*phip(:,2).*del2phi_L(:,2)-2/7*phip(:,3).*del2phi_L(:,3);
force_th_Lc(:,2)=-phi(:,2).*del2phi_L(:,1)+1/5*phi(:,2).*del2phi_L(:,3)+...
    3/5*phi(:,3).*del2phi_L(:,2);
force_th_Lc(:,3)=-phi(:,3).*del2phi_L(:,1)-1/3*phi(:,2).*del2phi_L(:,2)-...
    1/7*phi(:,3).*del2phi_L(:,3);
dr_force_th_Lc=zeros(N,n_modes);
dr_force_th_Lc(:,3)=D1*force_th_Lc(:,3);
dr_force_th_Lc(:,2)=-phip(:,2).*del2phi_L(:,1)-...
    phi(:,2).*drdel2phi_L(:,1).*(sign(r(floor(N/8))-r)+1)/2+...
    D1*(1/5*phi(:,2).*del2phi_L(:,3)+3/5*phi(:,3).*del2phi_L(:,2));
    %}

force_r_R=-gradphi_r_R.*del2phi_R;
force_th_R=-gradphi_th_R.*del2phi_R;
force_r_L=R2L(force_r_R,n_modes);
force_th_L=R2dL(force_th_R,n_modes);
r_force_th_L=diag(r)*force_th_L;
dr_force_th_L=D1*r_force_th_L;

rhs=force_r_L-dr_force_th_L;

bnf=1/2*gradphi_th_R(N,:).^2.*absc*wts'*3/2-force_th_L(N,2);
bss=R2dL(-gradphi_r_R(N,:).*gradphi_th_R(N,:),n_modes);

for j=n_modes-1:-1:1
    in1=IN\([0;rhs(2:N,j+1)].*r.^(j+3));
    in2=IN\([0;in1(2:N)]./r.^(2*j+4));
    in3=IN\([0;in2(2:N)].*r);
    in4p=in3.*r.^(2*j-2);
    in4=IN\([0;in4p(2:N)]);
    ps=in4.*r.^(2-j);
    if j>=2
        U2(j)=-(in2(N)+bss(j+1))/2/(2*j+1);
        U1(j)=-U2(j)-in4(N);
        sl(:,j+1)=ps+U2(j)./r.^j+U1(j)*r.^(2-j);
        slp(:,j+1)=(2-j)*r.^(1-j).*in4+r.^(2-j).*in4p-j*U2(j)./r.^(j+1)+(2-j)*U1(j)*r.^(1-j);
    else
        U2(j)=(in2(N)+bss(2))/6;
        %U3=(in1(N)+3*in2(N)-6*in3(N)-bnf-bss(2))/6;
        U3=U2(j)-in4(N);
        %U2=U1/30-in4(N);
        %U1(j)=U2(j)-in4(N)-U3;
        U1(j)=0;
        sl(:,2)=ps-U2(j)./r+U3*r.^2+U1(j)*r;
        slp(:,2)=in4+r.*in4p+U2(j)./r.^2+2*U3*2*r+U1(j);
        %sl(:,2)=ps-U1/30./r+U2*r.^2;
    end
end
[u_R,v_R]=Calc_u_v(r,sl,slp);
[cp_r_R,cp_th_R]=Calc_grad_R(r,cp2,D1);
[cm_r_R,cm_th_R]=Calc_grad_R(r,cm2,D1);
Lcp_R=u_R.*cp_r_R+v_R.*cp_th_R-par_NP*...
    (cp_r_R.*gradphi_r_R+cp_th_R.*gradphi_th_R+(cp2*P_mx(1:n_modes,:)).*del2phi_R);
Lcm_R=u_R.*cm_r_R+v_R.*cm_th_R+par_NP*...
    (cm_r_R.*gradphi_r_R+cm_th_R.*gradphi_th_R+(cm2*P_mx(1:n_modes,:)).*del2phi_R);
Lcp=R2L(Lcp_R,n_modes);
Lcm=R2L(Lcm_R,n_modes);
for j=1:n_modes
    Lcp(:,j)=Lcp(:,j)-par_NP*(D2+diag(2./r)*D1-diag(j*(j-1)./r.^2))*cp2(:,j);
    Lcm(:,j)=Lcm(:,j)-par_NP*(D2+diag(2./r)*D1-diag(j*(j-1)./r.^2))*cm2(:,j);
end
Lcp_L=Lcp;
Lcm_L=Lcm;

%Lcp(1,1)=cm2(1,1)-cp2(1,1);
%Lcm(1,1)=cm2(1,1)+cp2(1,1)-2;
Lcm(1,1)=cm2(1,1)-1;
Lcp(1,1)=cp2(1,1)-1;
Lcm(1,2)=D1(1,:)*cm2(:,2)+2/r(1)*cm2(1,2);
Lcp(1,2)=D1(1,:)*cp2(:,2)+2/r(1)*cp2(1,2);
Lcm(1,3)=D1(1,:)*cm2(:,3)+3/r(1)*cm2(1,3);
Lcp(1,3)=D1(1,:)*cp2(:,3)+3/r(1)*cp2(1,3);
%Lcm(1,3)=cm2(1,3);
%Lcp(1,3)=cp2(1,3);
Lcm(N,:)=D1(N,:)*cm2+par_CD/par_PE*cm2(N,:);
Lcp(N,:)=D1(N,:)*cp2-par_CD/par_PE*cp2(N,:);
%Lcm(N,:)=D1(N,:)*(cm2+cp2)-par_CD/par_PE*(cp2(N,:)-cm2(N,:));
%Lcp(N,:)=D1(N,:)*(cp2-cm2)-par_CD/par_PE*(cp2(N,:)+cm2(N,:));
%Lcm(N,1)=D1(N,:)*cm2(:,1)+par_CD/par_PE*cm2(N,1);
%Lcp(N,1)=D1(N,:)*cp2(:,1)-par_CD/par_PE*cp2(N,1);
Fres=reshape([Lcm;Lcp],[],1);
end

