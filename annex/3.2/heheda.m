clf

k_min=0.02;
k_max=0.18;
n_k=30;
s_min=2.8;
s_max=5;
n_s=30;

k = linspace(k_min, k_max, n_k);
s = linspace(s_min, s_max, n_s);
[K, S] = meshgrid(k, s);

a_k=0.004;
k_0=0.2;
k_1=0.222;
n=2;
p=5;
b_k=0.08;
b_s=0.8;
parameter=[a_k,k_0,k_1,n,p,b_k,b_s];

U = a_k+b_k.*K.^n./(k_0.^n+K.^n)-K./(1+K+S);
V = b_s./(1+(K./k_1).^p)-S./(1+K+S);
% quiver(K, S, U, V, 3,'ShowArrowHead','off','Marker','.','Color','#0072BD')
xlim([k_min k_max])
ylim([s_min s_max])
hold on

fun = @(x)paramfun(x,parameter);

x0_1 = [0.029,4.11];
x_1 = fsolve(fun,x0_1);

x0_2 = [0.08,4.177];
x_2 = fsolve(fun,x0_2);

x0_3 = [0.12,3.65];
x_3 = fsolve(fun,x0_3);

plot(x_1(1),x_1(2),'.b',x_2(1),x_2(2),'.g',x_3(1),x_3(2),'.r','MarkerSize',48)
hold on

% total=2000;
% k0=0.0931;
% s0=4.586;
% 
% tspan = [1 total];
% x0 = [k0 s0]; % initial
% opts = odeset('RelTol',1e-6,'AbsTol',1e-5);
% [t,y] = ode45(@(t,x) heheda_ode(t,x,parameter), tspan, x0,opts);
% plot(y(1:end,1),y(1:end,2))

k0_min=0.055;
k0_max=0.085;
n_k0=10;
s0_min=4.3;
s0_max=5;
n_s0=10;

ylist_exc=[];
ylist_non=[];
total=2000;
tspan=[1 total];
% opts = odeset('MaxStep', 1, 'RelTol',1e-2,'AbsTol',1e-4);
opts = odeset('MaxStep', 1);
for m=linspace(k0_min,k0_max,n_k0)
    for n=linspace(s0_min,s0_max,n_s0)
        x0=[m n];
        [t,y] = ode45(@(t,x) heheda_ode(t,x,parameter), tspan, x0,opts);
        if max(y(1:end,1))>0.12
            ylist_exc=[ylist_exc,y];
        else
            ylist_non=[ylist_non,y];
        end
    end
end

%%
for p=1:size(ylist_exc,2)/2
    plot(ylist_exc(:,2.*p-1),ylist_exc(:,2.*p),'Color','#D95319','linewidth',1.2)
    hold on
end

for q=1:size(ylist_non,2)/2
    plot(ylist_non(:,2.*q-1),ylist_non(:,2.*q),'Color','#4DBEEE','linewidth',1.2)
    hold on
end
%%
for i=1:size(ylist_exc,2)/2
    %plot(t,ylist_exc(:,2*i-1),'linewidth',1.2)
    xlim([0.02 0.18])
    ylim([2.8 5])
    xlabel('[ComK]','FontSize',20)
    ylabel('[ComS]','FontSize',20)
    set(gca,'FontSize',20)
    hold on
    text(0.025,4.3,'Node','FontSize',20)
    text(0.075,4.1,'Saddle','FontSize',20)
    text(0.11,3.5,'Unstable Focus','FontSize',20)
end