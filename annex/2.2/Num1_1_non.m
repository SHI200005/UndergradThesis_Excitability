% Num1 给出特征值实部
clf % 清除之前的图
tic

left=0.05;
bottom=0.18;
width=0.26;
height=0.72;
inter_w=0.06;
inter_h=0.12;

w=[0,1,2];
% h=[1,0];
h=0;

p_w=left+w.*(width+inter_w);
p_h=bottom+h.*(height+inter_h);

list=readmatrix('Num1_1_se.xlsx','Sheet','1','Range','A1:M6');
for ii=4 % for ii=1:2
    % Default Parameter
    %parameter=list(i,2:end);
    n=2;% 固定不变
    theta_y=list(ii,2);
    theta_z=list(ii,3);
    K=list(ii,4);
    % eta=0.01-100 参数范围
    kappa=list(ii,5);
    alpha=list(ii,6);
    beta_x=list(ii,7);
    beta_y=list(ii,8);
    beta_z=list(ii,9);
    
    Naccurate=1e3;%参数eta的精度
    fun_max=@(x)1e2.*(kappa.*x.^n./(1+x.^n)+beta_y).^n./(1+(kappa.*x.^n./(1+x.^n)+beta_y).^n+(alpha.*((kappa.*x.^n./(1+x.^n)+beta_y)./K).^n...,
        ./(1+((kappa.*x.^n./(1+x.^n)+beta_y)./K).^n)+beta_z).^n)+beta_x-x;
    x_max=fzero(fun_max,0.1);%初始值一定要取小
    x=linspace(beta_x,x_max,Naccurate);
    
    % Jacobian Matrix Element
    b=(x-beta_x).*n.*(1+(alpha.*((kappa.*x.^n./(1+x.^n)+beta_y)./K).^n./(1+((kappa.*x.^n./(1+x.^n)+beta_y)./K).^n)+beta_z).^n)...,
        ./(kappa.*x.^n./(1+x.^n)+beta_y)./(1+(kappa.*x.^n./(1+x.^n)+beta_y).^n+(alpha.*((kappa.*x.^n./(1+x.^n)+beta_y)./K).^n...,
        ./(1+((kappa.*x.^n./(1+x.^n)+beta_y)./K).^n)+beta_z).^n);
    c=-(x-beta_x).*n.*(alpha.*((kappa.*x.^n./(1+x.^n)+beta_y)./K).^n./(1+((kappa.*x.^n./(1+x.^n)+beta_y)./K).^n)+beta_z).^(n-1)...,
        ./(1+(kappa.*x.^n./(1+x.^n)+beta_y).^n+(alpha.*((kappa.*x.^n./(1+x.^n)+beta_y)./K).^n./(1+((kappa.*x.^n./(1+x.^n)+beta_y)./K).^n)+beta_z).^n);
    d=kappa.*n.*x.^(n-1)/theta_y./(1+x.^n).^2;
    h=alpha.*n.*((kappa.*x.^n./(1+x.^n)+beta_y)./K).^(n-1)./theta_z./K./(1+((kappa.*x.^n./(1+x.^n)+beta_y)./K).^n).^2;
    
    a=-ones(1,Naccurate);
    e=-ones(1,Naccurate)./theta_y;
    f=zeros(1,Naccurate);
    g=zeros(1,Naccurate);
    i=-ones(1,Naccurate)./theta_z;
    
    J=[a;b;c;d;e;f;g;h;i];
    
    % Calculate the eigan value
    A=zeros(Naccurate,3);
    for i=1:size(A,1)
        A(i,:)=eig([J(1,i),J(2,i),J(3,i);J(4,i),J(5,i),J(6,i);J(7,i),J(8,i),J(9,i)]);
    end
    
    eta=(x-beta_x).*(1+(kappa.*x.^n./(1+x.^n)+beta_y).^n+(alpha.*((kappa.*x.^n./(1+x.^n)+beta_y)./K).^n...,
        ./(1+((kappa.*x.^n./(1+x.^n)+beta_y)./K).^n)+beta_z).^n)./(kappa.*x.^n./(1+x.^n)+beta_y).^n;

   %%
    if ii==1
        
        subplot('Position',[p_w(1) p_h(1) width height])
        plot(eta,x,'k','linewidth',1.5)
        set(gca,'FontSize',12)
        title('1-Parameter Diagram','FontSize',14)
        xlabel('\eta','FontSize',14)
        ylabel('[x]','FontSize',14)
        
        subplot('Position',[p_w(2) p_h(1) width height])
        plot(x(2:end),real(A(2:end,1)),'r',x(2:end),real(A(2:end,2)),'g',x(2:end),real(A(2:end,3)),'b','linewidth',1.5)
        xlim([0 87])
        ylim([-1.2 0.2])
        set(gca,'YTick',-1.2:0.3:0.2,'FontSize',12)
        title('Eigenvalue','FontSize',14)
        xlabel('[x]','FontSize',14)
        ylabel('real(λ)','FontSize',14)
        
        [t,S,y]=time(2,1,80,30,50,n,theta_y,theta_z,K,kappa,alpha,beta_x,beta_y,beta_z);
        subplot('Position',[p_w(3) p_h(1) width height])
        yyaxis left
        plot(t,y(:,1),'linewidth',1.5)
        ylim([1.7 3])
        set(gca,'YTick',1.8:0.3:3,'FontSize',12)
        xlabel('Time','FontSize',14)
        ylabel('[x](\tau)','FontSize',14)
        yyaxis right
        plot(t,S,'linewidth',1)
        ylabel('\eta(\tau)','FontSize',14)
        ylim([1.9 3.1])
        set(gca,'YTick',1.9:0.3:3.1)
    end
    %%
    if ii==2
        IDsolution_SN(1)=list(ii,11);
        IDsolution_SN(2)=list(ii,12);
        
        subplot('Position',[p_w(1) p_h width height])
        plot(eta(1:IDsolution_SN(1)),x(1:IDsolution_SN(1)),'k',eta(IDsolution_SN(1):IDsolution_SN(2)),x(IDsolution_SN(1):IDsolution_SN(2)),'--k',eta(IDsolution_SN(2):end),x(IDsolution_SN(2):end),'k','linewidth',1.5)
        hold on
        plot(eta(IDsolution_SN(1)),x(IDsolution_SN(1)),'Color','b','Marker','.','MarkerSize',24)
        hold on
        plot(eta(IDsolution_SN(2)),x(IDsolution_SN(2)),'Color','b','Marker','.','MarkerSize',24)
        xlim([0 50])
        ylim([0 5.8])
        text(33,0.5,'SN 1')
        text(13,1.4,'SN 2')
        set(gca,'FontSize',12)
        title('1-Parameter Diagram','FontSize',14)
        xlabel('\eta','FontSize',14)
        ylabel('[x]','FontSize',14)        
        
        subplot('Position',[p_w(2) p_h width height])
        plot(x(2:end),real(A(2:end,2)),'g',x(2:end),real(A(2:end,3)),'b','linewidth',1.5)
        hold on
        plot(x(IDsolution_SN(1)),real(A(IDsolution_SN(1),2)),'Color','b','Marker','.','MarkerSize',24)
        hold on
        plot(x(IDsolution_SN(2)),real(A(IDsolution_SN(2),2)),'Color','b','Marker','.','MarkerSize',24)
        xlim([0,5.8])
        title('Eigenvalue')
        xlabel('[x]')
        ylabel('real(λ)')
        ylim([-1 0.5])
        text(0.15,-0.13,'SN 1')
        text(1.2,-0.13,'SN 2')
        set(gca,'FontSize',12)
        title('Eigenvalue','FontSize',14)
        xlabel('[x]','FontSize',14)
        ylabel('real(λ)','FontSize',14)        
        
        [t,S,y]=time(25,13,80,30,50,n,theta_y,theta_z,K,kappa,alpha,beta_x,beta_y,beta_z);
        subplot('Position',[p_w(3) p_h width height])
        yyaxis left
        plot(t,y(:,1),'linewidth',1.5)
        xlabel('Time','FontSize',14)
        ylabel('[x](t)','FontSize',14)
        xlim([0 80])
        ylim([0 4])
        xlabel('Time','FontSize',14)
        ylabel('[x](\tau)','FontSize',14)
        set(gca,'FontSize',12)
        yyaxis right
        plot(t,S,'linewidth',1)
        ylabel('\eta(\tau)','FontSize',14)     
        ylim([22 40])

    end
    %%
    if ii==3
        IDsolution_HB(1)=list(ii,10);
        IDsolution_HB(2)=list(ii,13);
    
        subplot('Position',[p_w(1) p_h(1) width height])
        plot(eta(1:IDsolution_HB(1)),x(1:IDsolution_HB(1)),'k',eta(IDsolution_HB(1):IDsolution_HB(2)),x(IDsolution_HB(1):IDsolution_HB(2)),'--k',eta(IDsolution_HB(2):end),x(IDsolution_HB(2):end),'k','linewidth',1.5)
        hold on
        plot(eta(IDsolution_HB(1)),x(IDsolution_HB(1)),'Color','r','Marker','.','MarkerSize',24)
        hold on
        plot(eta(IDsolution_HB(2)),x(IDsolution_HB(2)),'Color','r','Marker','.','MarkerSize',24)
        xlim([0 70])
        ylim([0 1.7])
        text(19,0.3,'HB 1')
        text(56,1,'HB 2')
        set(gca,'FontSize',12)
        title('1-Parameter Diagram','FontSize',14)
        xlabel('\eta','FontSize',14)
        ylabel('[x]','FontSize',14)
        
        subplot('Position',[p_w(2) p_h(1) width height])
        plot(x(2:end),real(A(2:end,2)),'g',x(2:end),real(A(2:end,3)),'b','linewidth',1.5)
        hold on
        plot(x(IDsolution_HB(1)),real(A(IDsolution_HB(1),2)),'Color','r','Marker','.','MarkerSize',24)
        hold on
        plot(x(IDsolution_HB(2)),real(A(IDsolution_HB(2),2)),'Color','r','Marker','.','MarkerSize',24)
        xlim([0 1.7])
        ylim([-1 0.5])
        text(0.3,-0.1,'HB 1')
        text(1,-0.1,'HB 2')  
        set(gca,'FontSize',12)
        title('Eigenvalue','FontSize',14)
        xlabel('[x]','FontSize',14)
        ylabel('real(λ)','FontSize',14)
        
        [t,S,y]=time(10,20,400,150,250,n,theta_y,theta_z,K,kappa,alpha,beta_x,beta_y,beta_z);
        subplot('Position',[p_w(3) p_h(1) width height])
        yyaxis left
        plot(t,y(:,1),'linewidth',1.5)
        xlabel('Time')
        ylabel('[x](t)')
        set(gca,'YTick',0:1:3,'FontSize',12)
        xlabel('Time','FontSize',14)
        ylabel('[x](\tau)','FontSize',14)
        yyaxis right
        plot(t,S,'linewidth',1)
        ylabel('\eta(\tau)','FontSize',14)
        ylim([7 33])
    end
    %%
    if ii==4
        IDsolution_HB(1)=list(ii,10);
        IDsolution_SN(1)=list(ii,11);
        IDsolution_SN(2)=list(ii,12);
        IDsolution_HB(2)=list(ii,13);
    
        subplot('Position',[p_w(1) p_h width height])
        plot(eta(1:IDsolution_HB(1)),x(1:IDsolution_HB(1)),'k',eta(IDsolution_HB(1):IDsolution_HB(2)),x(IDsolution_HB(1):IDsolution_HB(2)),'--k',eta(IDsolution_HB(2):end),x(IDsolution_HB(2):end),'k','linewidth',1.5)
        hold on
        plot(eta(IDsolution_HB(1)),x(IDsolution_HB(1)),'Color','r','Marker','.','MarkerSize',24)
        hold on
        plot(eta(IDsolution_SN(1)),x(IDsolution_SN(1)),'Color','b','Marker','.','MarkerSize',24)
        hold on
        plot(eta(IDsolution_SN(2)),x(IDsolution_SN(2)),'Color','b','Marker','.','MarkerSize',24)
        hold on
        plot(eta(IDsolution_HB(2)),x(IDsolution_HB(2)),'Color','r','Marker','.','MarkerSize',24)
        xlim([0 70])
        ylim([0 3.6])
        text(13,0.5,'HB 1')
        text(52,0.7,'SN 1')        
        text(34,1.3,'SN 2') 
        text(34,1.6,'HB 2')
        set(gca,'FontSize',12)
        xlabel('\eta','FontSize',14)
        ylabel('[x]','FontSize',14)
        title('1-Parameter Diagram','FontSize',14)
         
        subplot('Position',[p_w(2) p_h width height])
        plot(x(2:end),real(A(2:end,2)),'g',x(2:end),real(A(2:end,3)),'b','linewidth',1.5)
        hold on
        plot(x(IDsolution_HB(1)),real(A(IDsolution_HB(1),3)),'Color','r','Marker','.','MarkerSize',24)
        hold on
        plot(x(IDsolution_SN(1)),real(A(IDsolution_SN(1),3)),'Color','b','Marker','.','MarkerSize',24)
        hold on
        plot(x(IDsolution_SN(2)),real(A(IDsolution_SN(2),3)),'Color','b','Marker','.','MarkerSize',24)
        hold on
        plot(x(IDsolution_HB(2)),real(A(IDsolution_HB(2),3)),'Color','r','Marker','.','MarkerSize',24)
        xlim([0 3.6])
        ylim([-1 1])
        text(0,0.2,'HB 1')
        text(0.4,-0.15,'SN 1')        
        text(0.9,-0.15,'SN 2') 
        text(1.5,0.1,'HB 2')
        set(gca,'FontSize',12)
        title('Eigenvalue','FontSize',14)
        xlabel('[x]','FontSize',14)
        ylabel('real(λ)','FontSize',14)
        
        [t,S,y]=time(10,20,1600,600,1000,n,theta_y,theta_z,K,kappa,alpha,beta_x,beta_y,beta_z);
        subplot('Position',[p_w(3) p_h width height])
        yyaxis left
        plot(t,y(:,1),'linewidth',1.5)
        ylim([-1 14])
        set(gca,'ytick',0:4:14,'FontSize',12)
        xlabel('Time','FontSize',14)
        ylabel('[x](\tau)','FontSize',14)        
        yyaxis right
        plot(t,S,'linewidth',1)
        ylim([7 33])
        % set(gca,'FontSize',10)
        ylabel('\eta(\tau)','FontSize',14)
    end
end
toc

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