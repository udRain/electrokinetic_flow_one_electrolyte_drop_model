function data_r=tool_cpm_plot( r,cp,cm,n_th,label_opt,cur_time_i,par_dt )
%TOOL_PHI_PLOT Summary of this function goes here
%   Detailed explanation goes here
[N,M]=size(cp);
data=(cp+cm)/2;
tP_mx=zeros(M,16);
dth=2*pi/(n_th-1);
theta=0:dth:2*pi;
x=zeros(N,n_th);
y=zeros(N,n_th);
for j=1:n_th
    [tP_mx(:,j),~]...
        =Calc_LP(M-1,cos(theta(j)));
    for i=1:N
        x(i,j)=r(i)*cos(theta(j));
        y(i,j)=r(i)*sin(theta(j));
    end
end
data_r=data*tP_mx;
plot(x(N,:),y(N,:),'b');
hold on
%[C,h]=contour(x,y,data_r,[0.9991,0.9994,0.9997,0.9999,0.99999,1.00001,1.0001,1.0003,1.0006,1.0009],'b');
[C,h]=contour(x,y,data_r,[0.995,0.9985,0.9995,0.9999,1.0001,1.0005,1.0015,1.005],'b');
h.LineWidth=1;
%contour(x,y,data_r,1100,'b');
xlabel('x');
ylabel('y');
title(['t=' num2str(cur_time_i*par_dt)]);
if label_opt==1
    %clabel(C,h,'LabelSpacing',500);
    clabel(C,h,'manual','FontSize',20);
end
jpg_name=strcat('fig_t_',num2str(cur_time_i),'.fig');
saveas(gcf,jpg_name);
hold off
drawnow

%colormap jet;
%caxis([lr,ur]);
%view(0,0);
end

