clf
clear

list=readmatrix('Num1_1_se.xlsx','Sheet','2','Range','A1:R5');
ii=5;

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
x_s=0.27;
x_c=0.7;
eta_s=eta(x_s); % steady
eta_c=eta(x_c); % close/stimulus_1

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
A_s=eig(J_s);
A_c=eig(J_c);

[t_1,S_1,y_1]=time_1(eta_s,eta_c-eta_s,120,20,32,80,92,n,theta_y,theta_z,K,kappa,alpha,beta_x,beta_y,beta_z);

yyaxis left
plot(t_1,y_1(:,1),'linewidth',2.5)
xlabel('Time')
ylabel('[x](\tau)')
%ylim([0 10])
yyaxis right
plot(t_1,S_1,'linewidth',1.5)
ylabel('\eta(\tau)')
%ylim([40 90])
set(gca,'FontSize',16)
set(gca,'YTick',4.7:0.05:5.1)

%%
function[t,S,y]=time_1(eta0,S0,total,t_1,t_2,t_3,t_4,n,theta_y,theta_z,K,kappa,alpha,beta_x,beta_y,beta_z)

init_fun=@(x,eta)eta.*(kappa.*x.^n./(1+x.^n)+beta_y).^n./(1+(kappa.*x.^n./(1+x.^n)+beta_y).^n+(alpha.*((kappa.*x.^n./(1+x.^n)+beta_y)./K).^n./(1+((kappa.*x.^n./(1+x.^n)+beta_y)./K).^n)+beta_z).^n)+beta_x-x;

fun=@(x)init_fun(x,eta0);
x_init=fzero(fun,0);
y_init=kappa.*x_init.^n./(1+x_init.^n)+beta_y;
z_init=alpha*(y_init./K).^n./(1+(y_init./K).^n)+beta_z;

St=linspace(1,total,total);
S=S0.*(heaviside(St-t_1).*heaviside(t_4-St)-heaviside(St-t_2).*heaviside(t_3-St))+eta0;

tspan = [1 total];
y0 = [x_init y_init z_init];
opts = odeset('RelTol',1e-7,'AbsTol',1e-6);
[t,y] = ode45(@(t,y) Num1_1_ode(t,y,n,theta_y,theta_z,K,kappa,alpha,beta_x,beta_y,beta_z,St,S), tspan, y0,opts);

S=S0.*(heaviside(t-t_1).*heaviside(t_4-t)-heaviside(t-t_2).*heaviside(t_3-t))+eta0;

end


