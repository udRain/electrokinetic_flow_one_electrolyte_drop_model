clear all;
format long;
global P_mx;        %matrix of Legendre polynomial P_j
global dP_mx;       %matrix of 1st derivative of Legendre polynomial dP_j/dtheta
global absc;        %Abscissae of Gauss Legendre quadrature
global wts;         %Weights of Gauss Legendre quadrature
%{
v=VideoWriter('c_level_curve.mp4','MPEG-4');
v.Quality=100;
v.FrameRate=40;
open(v);
%}
setup_modes;
load('ic_lt.mat');
N=length(r);
par_PE=para(1);
par_NP=para(2);
par_CD=para(3);
par_EP=1;
par_dt=2e-7; 
max_t=1e-3;
iters=floor(max_t/par_dt);
ilevel=floor(iters/10);
t_v=0:par_dt:par_dt*iters;
n_modes=3;

I=eye(N);
[NPm_mx,NPp_mx]=tNP_lhs(r,D1,D2,n_modes,par_dt,para);
IN=D1;
IN(1,:)=I(1,:);

phi=zeros(N,n_modes);
cp=zeros(N,n_modes);
cm=zeros(N,n_modes);
sl=zeros(N,n_modes);
slp=zeros(N,n_modes);

phi(:,1)=phi_0;
cp(:,1)=cp_0;
cm(:,1)=cm_0;

dphi=zeros(N,n_modes);
d2phi_L=1/2/par_PE^2*(cm-cp);
d2phi_R=d2phi_L*P_mx(1:n_modes,:);
phi(:,2)=phi(:,2)-par_EP*(r)-par_EP/2./r.^2;
[gphi_r_R,gphi_th_R]=Calc_grad_R(r,phi,D1);
u_R=zeros(N,16);
v_R=zeros(N,16);

c1_b=zeros(1,n_modes);
U1=zeros(1,n_modes-1);
U2=zeros(1,n_modes-1);

Ufar=zeros(1,iters+1);
crs=zeros(1,iters+1);
U1rec=zeros(1,iters+1);

for i=1:iters
    
    [rhscm,rhscp]=tNP_rhs(r,cm,cp,gphi_r_R,gphi_th_R,d2phi_R,u_R,v_R,par_dt,par_NP,D1 );
    %for j=1:floor(2*M/3)
    for j=1:n_modes
        %Ncmj_mx=NPm_mx(:,:,j)+0*par_NP*(diag(dphitdr(:,1))*D1+diag(del2phi_L(:,1)));
        %Ncpj_mx=NPp_mx(:,:,j)-0*par_NP*(diag(dphitdr(:,1))*D1+diag(del2phi_L(:,1)));
        %Ncmj_mx(N,:)=D1(N,:)+par_CD/par_PE*I(N,:);
        %Ncmj_mx(1,:)=I(1,:);
        %Ncpj_mx(N,:)=D1(N,:)-par_CD/par_PE*I(N,:);
        %Ncpj_mx(1,:)=I(1,:);
        cm(:,j)=NPm_mx(:,:,j)\rhscm(:,j);
        cp(:,j)=NPp_mx(:,:,j)\rhscp(:,j);
    end
    %{
    plot(r,cm(:,3));
    hold on
    plot(r,cp(:,3));
    hold off
    drawnow
    %}
    d2phi_L=1/2*(cm-cp)/par_PE^2;
    d2phi_R=d2phi_L*P_mx(1:n_modes,:);
    rhspeb=zeros(1,n_modes);
    rhspeb(1)=-par_CD/par_PE;
    rhspeb(2)=par_EP;

    for j=1:n_modes
        in1=IN\(r.^(j+1).*[0;d2phi_L(2:N,j)]);
        in2p=r.^(-2*j).*in1;
        in2=IN\([0;in2p(2:N)]);
        c1_b(j)=(rhspeb(j)-in1(N)-(j-1)*in2(N))/j;
        phi(:,j)=r.^(j-1).*in2-c1_b(j)./r.^j;
        dphi(:,j)=(j-1)*r.^(j-2).*in2+r.^(j-1).*in2p+j*c1_b(j)./r.^(j+1);
    end
    phi(:,2)=phi(:,2)-par_EP*(r);
    dphi(:,2)=dphi(:,2)-par_EP*ones(size(r));
    gphi_r_R=dphi*P_mx(1:n_modes,:);
    gphi_th_R=diag(1./r)*phi*dP_mx(1:n_modes,:);
    crs(i+1)=c1_b(3);
    %{
    plot(r,chebdifft(cm(:,1),1)/ra.*dphi(:,3)+...
        chebdifft(cm(:,2),1)/ra.*dphi(:,2)/3*2*0+...
        chebdifft(cm(:,3),1)/ra.*dphi(:,1)+...
        chebdifft(cm(:,3),1)/ra.*dphi(:,3)/7*2*0);
    hold on
    plot(r,chebdifft(cp(:,1),1)/ra.*dphi(:,3)+...
        chebdifft(cp(:,2),1)/ra.*dphi(:,2)/3*2*0+...
        chebdifft(cp(:,3),1)/ra.*dphi(:,1)+...
        chebdifft(cp(:,3),1)/ra.*dphi(:,3)/7*2*0);
    hold off
    drawnow;
    %}
    
    force_r_R=-gphi_r_R.*d2phi_R;
    force_th_R=-gphi_th_R.*d2phi_R;
    force_r_L=R2L(force_r_R,n_modes);
    force_th_L=R2dL(force_th_R,n_modes);
    r_force_th_L=diag(r)*force_th_L;
    dr_force_th_L=D1*r_force_th_L;
    rhs=force_r_L-dr_force_th_L;
   
    bnf=1/2*gphi_th_R(N,:).^2.*absc*wts'*3/2-force_th_L(N,2);

    bss=R2dL(-gphi_r_R(N,:).*gphi_th_R(N,:),n_modes);

    for j=n_modes-1:-1:1
        in1s=IN\([0;rhs(2:N,j+1)].*r.^(j+3));
        in2s=IN\([0;in1s(2:N)]./r.^(2*j+4));
        in3s=IN\([0;in2s(2:N)].*r);
        in4p=in3s.*r.^(2*j-2);
        in4s=IN\([0;in4p(2:N)]);
        ps=in4s.*r.^(2-j);
        if j>=2
            U2(j)=-(in2s(N)+bss(j+1))/2/(2*j+1);
            U1(j)=-U2(j)-in4s(N);
            sl(:,j+1)=ps+U2(j)./r.^j+U1(j)*r.^(2-j);
            slp(:,j+1)=(2-j)*r.^(1-j).*in4s+r.^(2-j).*in4p-j*U2(j)./r.^(j+1)+(2-j)*U1(j)*r.^(1-j);
        else
            U2(j)=(in2s(N)+bss(2))/6;
            U3=(in1s(N)+3*in2s(N)-6*in3s(N)-bnf-bss(2))/6;
            Ufar(i+1)=U3*2;
            U1(j)=U2(j)-in4s(N)-U3;
            sl(:,2)=ps-U2(j)./r+U3*r.^2+U1(j)*r;
            slp(:,2)=in4s+r.*in4p+U2(j)./r.^2+2*U3*2*r+U1(j);
        end
    end
    U1rec(i+1)=U1(2);

    [u_R,v_R,~,~]=Calc_u_v(r,sl,slp);

    if i==floor(iters/4)
        display('25% completed.');
    else if i==floor(iters/2)
            display('50% completed.');
        else if i==floor(iters/4)*3
                display('75% completed.')
            end
        end
    end
    
    if mod(i,ilevel)==0
        data_r=tool_cpm_plot(r,cp,cm,128,i/ilevel>0,i,par_dt);
        %imv=getframe(gcf);
    end
end
%close(v);
para=[par_PE,par_NP,par_CD,par_EP];
%save('ic_0s.mat','r','ra','para','phi','cm','cp','sl','Ufar','-v7.3');