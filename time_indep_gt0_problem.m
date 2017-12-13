clear all;
format long;
global P_mx;        %matrix of Legendre polynomial P_j
global dP_mx;       %matrix of 1st derivative of Legendre polynomial dP_j/dtheta
global absc;        %Abscissae of Gauss Legendre quadrature
global wts;         %Weights of Gauss Legendre quadrature

setup_modes;
n_modes=3;
load('ic_lt_ss.mat');
N=length(r);
par_PE=para(1);
par_NP=para(2);
par_CD=para(3);
par_zeta=para(4);
par_EP=0.01;


I=eye(N);
IN=D1;
IN(1,:)=I(1,:);

eps=1e-7;

cm_0s=zeros(N,n_modes);
cp_0s=zeros(N,n_modes);
cm_0s(:,1)=cm_0;
cp_0s(:,1)=cp_0;
cm_0=cm_0s;
cp_0=cp_0s;


for i=1:100
    DF=zeros(2*N*n_modes);
    [Fu,~,~,~,~,phi_test]=F_cpm_g(r,n_modes,cm_0,cp_0,cm_0,cp_0,par_PE,par_CD,par_EP,par_NP,D1,D2);

    for j=1:N
        for k=1:n_modes
            cm_1=cm_0;
            cp_1=cp_0;
            cm_1(j,k)=cm_1(j,k)+eps;
            cp_1(j,k)=cp_1(j,k)+eps;
            tempv=(F_cpm_g(r,n_modes,cm_1,cp_0,cm_1,cp_0,par_PE,par_CD,par_EP,par_NP,D1,D2)-Fu)/eps;
            DF(:,(k-1)*2*N+j)=tempv;
            tempv=(F_cpm_g(r,n_modes,cm_0,cp_1,cm_0,cp_1,par_PE,par_CD,par_EP,par_NP,D1,D2)-Fu)/eps;
            DF(:,(k-1)*2*N+j+N)=tempv;
        end
    end
    du=-DF\Fu;
    un=[cm_0;cp_0]+1/2*reshape(du,2*N,n_modes);
    cm_0=un(1:N,:);
    cp_0=un(N+1:2*N,:);
    ea=norm(Fu)
    if ea<eps/par_PE^2/10
        break;
    end

end

[Fu,Lcp_R,Lcm_R,Lcp_L,Lcm_L]=F_cpm_g(r,n_modes,cm_0,cp_0,cm_0,cp_0,par_PE,par_CD,par_EP,par_NP,D1,D2);

rhspe=1/2*(cm_0-cp_0);
%rhspe(1:floor(N/5),1)=rhspe(1:floor(N/5),1)*0;
del2phi_L=rhspe/par_PE^2;

del2phi_R=del2phi_L*P_mx(1:n_modes,:);
rhspeb=zeros(1,n_modes);
rhspeb(1)=-par_CD/par_PE;
rhspeb(2)=par_EP;
phi=zeros(N,n_modes);
phip=zeros(N,n_modes);
sl=zeros(N,n_modes);
slp=zeros(N,n_modes);
c1_b=zeros(1,n_modes);
U1=zeros(1,n_modes-1);
U2=zeros(1,n_modes-1);

for j=1:n_modes
    inp1=IN\(r.^(j+1).*[0;del2phi_L(2:N,j)]);
    inp2p=r.^(-2*j).*inp1;
    inp2=IN\([0;inp2p(2:N)]);
    c1_b(j)=(rhspeb(j)-inp1(N)-(j-1)*inp2(N))/j;
    if j==1
        phi(:,j)=inp2;
        phip(:,j)=inp2p;
    else if j==2
            phi(:,j)=r.*inp2-par_EP/2./r.^2;
            %phi(:,j)=r.^(j-1).*inp2-c1_b./r.^j;
            phip(:,j)=(j-1)*r.^(j-2).*inp2+r.^(j-1).*inp2p+j*c1_b(j)./r.^(j+1);
        else
            phi(:,j)=r.^(j-1).*inp2-c1_b(j)./r.^j;
            phip(:,j)=(j-1)*r.^(j-2).*inp2+r.^(j-1).*inp2p+j*c1_b(j)./r.^(j+1);
        end
    end
    %phi(:,j)=r.^(j-1).*inp2-c1_b(j)./r.^j;
    %phip(:,j)=(j-1)*r.^(j-2).*inp2+r.^(j-1).*inp2p+j*c1_b(j)./r.^(j+1);
end
phi(:,2)=phi(:,2)-par_EP*r;
phip(:,2)=phip(:,2)-par_EP;
gradphi_r_R=phip*P_mx(1:n_modes,:);
gradphi_th_R=diag(1./r)*phi*dP_mx(1:n_modes,:);
[gradc_r_R,gradc_th_R]=Calc_grad_R(r,cm_0,D1);

force_r_Lc=zeros(N,n_modes);
force_th_Lc=zeros(N,n_modes);
%force_r_Lc(:,1)=-phip(:,1).*del2phi_L(:,1)-1/3*phip(:,2).*del2phi_L(:,2)-...
%    1/5*phip(:,3).*del2phi_L(:,3);
force_r_Lc(:,2)=-(phip(:,1).*del2phi_L(:,2)+phip(:,2).*del2phi_L(:,1))-...
    2/5*(phip(:,3).*del2phi_L(:,2)+phip(:,2).*del2phi_L(:,3));
force_r_Lc(:,3)=-(phip(:,1).*del2phi_L(:,3)+phip(:,3).*del2phi_L(:,1))-...
    2/3*phip(:,2).*del2phi_L(:,2)-2/7*phip(:,3).*del2phi_L(:,3);
force_th_Lc(:,2)=-phi(:,2).*del2phi_L(:,1)+1/5*phi(:,2).*del2phi_L(:,3)+...
    3/5*phi(:,3).*del2phi_L(:,2);
force_th_Lc(:,3)=-phi(:,3).*del2phi_L(:,1)-1/3*phi(:,2).*del2phi_L(:,2)-...
    1/7*phi(:,3).*del2phi_L(:,3);
dr_force_th_Lc=D1*force_th_Lc;

%{
force_r_R=-gradphi_r_R.*del2phi_R;
force_th_R=-gradphi_th_R.*del2phi_R;
force_r_L=R2L(force_r_R,n_modes);
force_th_L=R2dL(force_th_R,n_modes);
r_force_th_L=diag(r)*force_th_L;
dr_force_th_L=D1*r_force_th_L;
%}
rhs=force_r_Lc-dr_force_th_Lc;

%rhs=force_r_L-dr_force_th_L;
%rhs(1:floor(N/5) ,:)=rhs(1:floor(N/5),:)*0;

bnf=1/2*gradphi_th_R(N,:).^2.*absc*wts'*3/2-force_th_Lc(N,2);
bss=R2dL(-gradphi_r_R(N,:).*gradphi_th_R(N,:),n_modes);

for j=1:n_modes-1
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
        Ufar=U3*2;
        U1(j)=U2(j)-in4(N)-U3;
        sl(:,2)=ps-U2(j)./r+U3*r.^2+U1(j)*r;
        slp(:,2)=in4+r.*in4p+U2(j)./r.^2+2*U3*2*r+U1(j);
        %sl(:,2)=ps-U1/30./r+U2*r.^2;
    end
end
[u_R,v_R]=Calc_u_v(r,sl,slp);
save('ic_0s_ss.mat','r','ra','para','phi','cm_0','cp_0','sl','Ufar','-v7.3');