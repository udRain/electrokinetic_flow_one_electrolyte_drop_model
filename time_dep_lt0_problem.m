format long;
N=160;
n_rpts=N;
rmin=1;
rmax=2;
ra=(rmax-rmin)/2;
[x,DM]=chebdif(N,3);
D1=DM(:,:,1)/ra;
D2=DM(:,:,2)/ra^2;
D3=DM(:,:,3)/ra^3;
r=rmax-ra+ra*x;

dt=1e-7;
MAXIT=100000;

par_PE=1/300;
par_NP=5;
par_CD=2;
par_zeta=2*asinh(par_CD/2);
k=sqrt(1/par_NP/dt);

I=eye(N);
IN_mx=D1;
IN_mx(1,:)=I(1,:);

sc=exp(-(r-1)/par_PE);
phi_0=atanh(sc*tanh(par_zeta/4))*4;
cm_0=exp(phi_0);
cp_0=exp(-phi_0);

cm_0p=chebdifft(cm_0,1)/ra;
cp_0p=chebdifft(cp_0,1)/ra;

phi_s=phi_0;
phi_old=phi_0;
cm_s=cm_0;
cp_s=cp_0;

c1b=zeros(MAXIT,1);
dphidrb=-par_CD/par_PE;

for i=1:MAXIT
    phi_0pp=(cm_0-cp_0)/2/par_PE^2;
    
    in1=IN_mx\(r.^2.*[0;phi_0pp(2:N)]);
    c1_b=dphidrb-in1(N);
    c1b(i)=c1_b;
    in2=IN_mx\(1./r.^2.*[0;in1(2:N)]);
    phi_0=in2-0*c1_b./r;
    phi_0p=in1./r.^2;
    
    Ncm_mx=I/dt/par_NP-(D2+diag(2./r)*D1)+(diag(phi_0pp)+diag(phi_0p)*D1);
    Ncp_mx=I/dt/par_NP-(D2+diag(2./r)*D1)-(diag(phi_0pp)+diag(phi_0p)*D1);
    %Ncm_mx=I/dt/par_NP-(D2+diag(2./r)*D1);
    %Ncp_mx=I/dt/par_NP-(D2+diag(2./r)*D1);
    Ncm_mx(N,:)=D1(N,:)-dphidrb*I(N,:);
    %Ncm_mx(N,:)=D1(N,:);
    Ncm_mx(1,:)=I(1,:);
    Ncp_mx(N,:)=D1(N,:)+dphidrb*I(N,:);
    %Ncp_mx(N,:)=D1(N,:);
    Ncp_mx(1,:)=I(1,:);
    
    %rhs_cc=cc_0/dt/par_NP+cq_0p.*phi_0p+cq_0.*phi_0pp;
    %rhs_cq=cq_0/dt/par_NP+cc_0p.*phi_0p+cc_0.*phi_0pp;
    rhs_cm=cm_0/dt/par_NP;
    rhs_cp=cp_0/dt/par_NP;
    rhs_cm(1)=1;
    rhs_cp(1)=1;
    rhs_cm(N)=0;
    rhs_cp(N)=0;
    cm_0=Ncm_mx\rhs_cm;
    cp_0=Ncp_mx\rhs_cp;

    dif=norm(phi_0-phi_old)
    if dif<5e-7
        break;
    end
    phi_old=phi_0;
    %plot(r,cm_0-exp(phi_0));
    %drawnow;
    %{
    rhs_cmc=-cm_0c/dt/par_NP+cm_0p.*phi_0p+(cm_0c+1).*phi_0pp;
    rhs_cpc=-cp_0c/dt/par_NP-cp_0p.*phi_0p-(cp_0c+1).*phi_0pp;
    cmb=-dphidrb;
    cpb=dphidrb;
    cmin1=IN_mx\([0;rhs_cmc(2:N)].*r.*exp(k*(r-1)));
    cpin1=IN_mx\([0;rhs_cpc(2:N)].*r.*exp(k*(r-1)));
    %{
    for id=1:N
        if cmin1(id)<1e-12
            cmin1(id)=0;
        end
        if cpin1(id)<1e-12
            cpin1(id)=0;
        end
    end
    %}
    cmin2=IN_mx\(cmin1.*exp(-2*k*(r-1)));
    cpin2=IN_mx\(cpin1.*exp(-2*k*(r-1)));
    cmc=((1-k-cmb)*cmin2(N)-cmin1(N))/(1+k-cmb);
    cpc=((1-k-cpb)*cpin2(N)-cpin1(N))/(1+k-cpb);
    cm_0c=cmin2./r.*exp(k*(r-1))-cmc./r.*exp(-k*(r-1));
    cp_0c=cpin2./r.*exp(k*(r-1))-cpc./r.*exp(-k*(r-1));
    cm_0p=chebdifft(cm_0c,1)/ra;
    cp_0p=chebdifft(cp_0c,1)/ra;
    %}
end
%phi_0i=ra*(ws.*phi_0pp')*G_lp;
%phi_0=phi_0i'+par_CD/par_PE./r;
    
para=[par_PE,par_NP,par_CD];

save('ic_lt.mat','x','r','ra','para','phi_0','cm_0','cp_0','D1','D2','-v7.3');