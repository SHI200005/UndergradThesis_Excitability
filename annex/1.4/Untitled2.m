% Num1 给出特征值实部
clf
tic

% Default Parameter
%parameter=list(i,2:end);
n=2;% 固定不变
theta_y=1.010;
theta_z=71.57;
K=1.381;
% eta=0.01-100 参数范围
kappa=1.234;
alpha=12.95;
beta_x=0.04663;
beta_y=0.04206;
beta_z=0.2087;

Naccurate=1e3;%参数eta的精度
x_max=4.059;%初始值一定要取小
x=linspace(beta_x,x_max,Naccurate);

eta=(x-beta_x).*(1+(kappa.*x.^n./(1+x.^n)+beta_y).^n+(alpha.*((kappa.*x.^n./(1+x.^n)+beta_y)./K).^n...,
    ./(1+((kappa.*x.^n./(1+x.^n)+beta_y)./K).^n)+beta_z).^n)./(kappa.*x.^n./(1+x.^n)+beta_y).^n;

%%
IDsolution_SN(1)=25;
IDsolution_SN(2)=119;
IDsolution_HB(1)=381;

plot(eta(1:IDsolution_SN(1)),x(1:IDsolution_SN(1)),'k',eta(IDsolution_SN(1):IDsolution_HB(1)),x(IDsolution_SN(1):IDsolution_HB(1)),'--k',...,
    eta(IDsolution_HB(1):end),x(IDsolution_HB(1):end),'k','linewidth',1.5)
hold on
plot(eta(IDsolution_SN(1)),x(IDsolution_SN(1)),'Color','b','Marker','.','MarkerSize',24)
hold on
plot(eta(IDsolution_SN(2)),x(IDsolution_SN(2)),'Color','b','Marker','.','MarkerSize',24)
hold on
plot(eta(IDsolution_HB(1)),x(IDsolution_HB(1)),'Color','r','Marker','.','MarkerSize',24)
title('1-Parameter Diagram')
xlabel('\eta')
ylabel('[x]')
xlim([0 50])
ylim([0 2])
set(gca,'XTick',0:10:50)
set(gca,'YTick',0:0.4:2)
set(gca,'FontSize',12)

toc