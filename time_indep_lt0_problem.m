format long;
N=100;
n_rpts=N;
rmin=1;
rmax=1.2;
ra=(rmax-rmin)/2;
[x,DM]=chebdif(N,3);
D1=DM(:,:,1)/ra;
D2=DM(:,:,2)/ra^2;
r=rmax-ra+ra*x;

IN_mx=D1;
IN_mx(1,:)=IN_mx(1,:)*0;
IN_mx(1,1)=1;

par_PE=1/400;
par_NP=5;
par_CD=1;
par_zeta=2*asinh(par_CD/2);

sc=exp(-(r-1)/par_PE);
phi_s=atanh(sc*tanh(par_zeta/4))*4;
cm_0=exp(phi_s);
cp_0=exp(-phi_s);

%load('ic_lt.mat','cm_0','cp_0');
dphidrb=-par_CD/par_PE;
lbc=dphidrb;
eps=1e-9;

for i=1:1000
    DF=zeros(2*N);
    [Fu,phi]=F_cpm(r,cm_0,cp_0,lbc,par_PE,D1,D2);
    
    
    for j=1:N
        cm_1=cm_0;
        cp_1=cp_0;
        cm_1(j)=cm_1(j)+eps;
        cp_1(j)=cp_1(j)+eps;
        tempv=(F_cpm(r,cm_1,cp_0,lbc,par_PE,D1,D2)-Fu)/eps;
        DF(:,j)=tempv;
        tempv=(F_cpm(r,cm_0,cp_1,lbc,par_PE,D1,D2)-Fu)/eps;
        DF(:,j+N)=tempv;
    end
    du=-DF\Fu;
    un=[cm_0;cp_0]+du;
    cm_0=un(1:N);
    cp_0=un(N+1:2*N);
    dif=norm(Fu)

    if dif<1e-7
        break;
    end
end

para=[par_PE,par_NP,par_CD,par_zeta];
save('ic_lt_ss.mat','r','ra','para','phi','cm_0','cp_0','D1','D2','-v7.3');