function dydt=Num1_1_ode(t,y,n,theta_y,theta_z,K,kappa,alpha,beta_x,beta_y,beta_z,eta,St,S)
S=interp1(St,S,t);

dydt=zeros(3,1);
dydt(1)=S.*y(2).^n./(1+y(2).^n+y(3).^n)-y(1)+beta_x; % x
dydt(2)=(kappa.*y(1).^n/(1+y(1).^n)-y(2)+beta_y)./theta_y; % y
dydt(3)=(alpha*(y(2)./K).^n/(1+(y(2)./K).^n)-y(3)+beta_z)./theta_z; % z