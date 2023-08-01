function dxdt=heheda_ode(t,x,parameter)

a_k=parameter(1);
k_0=parameter(2);
k_1=parameter(3);
n=parameter(4);
p=parameter(5);
b_k=parameter(6);
b_s=parameter(7);

dxdt=zeros(2,1);
dxdt(1) = a_k+b_k.*x(1).^n./(k_0.^n+x(1).^n)-x(1)./(1+x(1)+x(2));
dxdt(2) = b_s./(1+(x(1)./k_1).^p)-x(2)./(1+x(1)+x(2));
