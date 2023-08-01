clf

left=0.04;
bottom=0.1;
width=0.25;
height=0.36;
inter_w=0.08;
inter_h=0.12;
w=[0,1,2];
h=[1,0];
p_w=left+w.*(width+inter_w);
p_h=bottom+h.*(height+inter_h);

list=readmatrix('Num1_1_se.xlsx','Sheet','2','Range','A1:R5');
ii=2;

n=2;% ¹Ì¶¨²»±ä
theta_y=list(ii,2);
theta_z=list(ii,3);
K=list(ii,4);
kappa=list(ii,5);
alpha=list(ii,6);
beta_x=list(ii,7);
beta_y=list(ii,8);
beta_z=list(ii,9);

eta=@(x)(x-beta_x).*(1+(kappa.*x.^n./(1+x.^n)+beta_y).^n+(alpha.*((kappa.*x.^n./(1+x.^n)+beta_y)./K).^n...,
    ./(1+((kappa.*x.^n./(1+x.^n)+beta_y)./K).^n)+beta_z).^n)./(kappa.*x.^n./(1+x.^n)+beta_y).^n;
x_s=0.085;
x_c=1.1;
x_f=1.18;
eta_s=eta(x_s); % steady
eta_c=eta(x_c); % close/stimulus_1
eta_f=eta(x_f); % far/stimulus_2

b=@(x)(x-beta_x).*n.*(1+(alpha.*((kappa.*x.^n./(1+x.^n)+beta_y)./K).^n./(1+((kappa.*x.^n./(1+x.^n)+beta_y)./K).^n)+beta_z).^n)...,
    ./(kappa.*x.^n./(1+x.^n)+beta_y)./(1+(kappa.*x.^n./(1+x.^n)+beta_y).^n+(alpha.*((kappa.*x.^n./(1+x.^n)+beta_y)./K).^n...,
    ./(1+((kappa.*x.^n./(1+x.^n)+beta_y)./K).^n)+beta_z).^n);
c=@(x)-(x-beta_x).*n.*(alpha.*((kappa.*x.^n./(1+x.^n)+beta_y)./K).^n./(1+((kappa.*x.^n./(1+x.^n)+beta_y)./K).^n)+beta_z).^(n-1)...,
    ./(1+(kappa.*x.^n./(1+x.^n)+beta_y).^n+(alpha.*((kappa.*x.^n./(1+x.^n)+beta_y)./K).^n./(1+((kappa.*x.^n./(1+x.^n)+beta_y)./K).^n)+beta_z).^n);
d=@(x)kappa.*n.*x.^(n-1)/theta_y./(1+x.^n).^2;
h=@(x)alpha.*n.*((kappa.*x.^n./(1+x.^n)+beta_y)./K).^(n-1)./theta_z./K./(1+((kappa.*x.^n./(1+x.^n)+beta_y)./K).^n).^2;
a=-1;
e=-1./theta_y;
f=0;
g=0;
i=-1./theta_z;

J_s=[a,b(x_s),c(x_s);d(x_s),e,f;g,h(x_s),i];
J_c=[a,b(x_c),c(x_c);d(x_c),e,f;g,h(x_c),i];
J_f=[a,b(x_f),c(x_f);d(x_f),e,f;g,h(x_f),i];
A_s=eig(J_s);
A_c=eig(J_c);
A_f=eig(J_f);

[t_1,S_1,y_1]=time(eta_s,eta_c-eta_s,60,20,28.1,n,theta_y,theta_z,K,kappa,alpha,beta_x,beta_y,beta_z);
[t_2,S_2,y_2]=time(eta_s,eta_c-eta_s,60,20,32,n,theta_y,theta_z,K,kappa,alpha,beta_x,beta_y,beta_z);
[t_3,S_3,y_3]=time(eta_s,eta_c-eta_s,250,30,215,n,theta_y,theta_z,K,kappa,alpha,beta_x,beta_y,beta_z);
[t_4,S_4,y_4]=time(eta_s,eta_f-eta_s,60,20,24,n,theta_y,theta_z,K,kappa,alpha,beta_x,beta_y,beta_z);
[t_5,S_5,y_5]=time(eta_s,eta_f-eta_s,60,20,26,n,theta_y,theta_z,K,kappa,alpha,beta_x,beta_y,beta_z);
[t_6,S_6,y_6]=time(eta_s,eta_f-eta_s,250,30,215,n,theta_y,theta_z,K,kappa,alpha,beta_x,beta_y,beta_z);

subplot('Position',[p_w(1) p_h(1) width height])
yyaxis left
plot(t_1,y_1(:,1),'linewidth',1.5)
xlabel('Time')
ylabel('[x](\tau)')
ylim([-1 6])
yyaxis right
plot(t_1,S_1,'linewidth',1)
ylabel('\eta(\tau)')
ylim([4 17])
set(gca,'FontSize',16)

subplot('Position',[p_w(2) p_h(1) width height])
yyaxis left
plot(t_2,y_2(:,1),'linewidth',1.5)
xlabel('Time')
ylabel('[x](\tau)')
ylim([-1 6])
yyaxis right
plot(t_2,S_2,'linewidth',1)
ylabel('\eta(\tau)')
ylim([4 17])
set(gca,'FontSize',16)

subplot('Position',[p_w(3) p_h(1) width height])
yyaxis left
plot(t_3,y_3(:,1),'linewidth',1.5)
xlabel('Time')
ylabel('[x](\tau)')
ylim([-1 6])
yyaxis right
plot(t_3,S_3,'linewidth',1)
ylabel('\eta(\tau)')
ylim([4 17])
xlim([0 250])
set(gca,'FontSize',16)

subplot('Position',[p_w(1) p_h(2) width height])
yyaxis left
plot(t_4,y_4(:,1),'linewidth',1.5)
xlabel('Time')
ylabel('[x](\tau)')
ylim([-1 6])
yyaxis right
plot(t_4,S_4,'linewidth',1)
ylabel('\eta(\tau)')
ylim([4 17])
set(gca,'FontSize',16)

subplot('Position',[p_w(2) p_h(2) width height])
yyaxis left
plot(t_5,y_5(:,1),'linewidth',1.5)
xlabel('Time')
ylabel('[x](\tau)')
ylim([-1 6])
yyaxis right
plot(t_5,S_5,'linewidth',1)
ylabel('\eta(\tau)')
ylim([4 17])
set(gca,'FontSize',16)

subplot('Position',[p_w(3) p_h(2) width height])
yyaxis left
plot(t_6,y_6(:,1),'linewidth',1.5)
xlabel('Time')
ylabel('[x](\tau)')
ylim([-1 6])
yyaxis right
plot(t_6,S_6,'linewidth',1)
ylabel('\eta(\tau)')
ylim([4 17])
xlim([0 250])
set(gca,'FontSize',16)

function[t,S,y]=time(eta0,S0,total,start,endtime,n,theta_y,theta_z,K,kappa,alpha,beta_x,beta_y,beta_z)

init_fun=@(x,eta)eta.*(kappa.*x.^n./(1+x.^n)+beta_y).^n./(1+(kappa.*x.^n./(1+x.^n)+beta_y).^n+(alpha.*((kappa.*x.^n./(1+x.^n)+beta_y)./K).^n./(1+((kappa.*x.^n./(1+x.^n)+beta_y)./K).^n)+beta_z).^n)+beta_x-x;

fun=@(x)init_fun(x,eta0);
x_init=fzero(fun,0);
y_init=kappa.*x_init.^n./(1+x_init.^n)+beta_y;
z_init=alpha*(y_init./K).^n./(1+(y_init./K).^n)+beta_z;

St=linspace(1,total,total);
S=S0.*heaviside(St-start).*heaviside(endtime-St)+eta0;

tspan = [1 total];
y0 = [x_init y_init z_init];
opts = odeset('RelTol',1e-6,'AbsTol',1e-5);
[t,y] = ode45(@(t,y) Num1_1_ode(t,y,n,theta_y,theta_z,K,kappa,alpha,beta_x,beta_y,beta_z,eta0,St,S), tspan, y0,opts);

S=S0.*heaviside(t-start).*heaviside(endtime-t)+eta0;

end


