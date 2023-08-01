clear
clf
tic

init_m_X=0;
init_m_Y=0;
init_m_Z=0;
init_X=290;
init_Y=570;
init_Z=870;
t=0;

% 参数
n=2;

K_yx=4e3;
K_zx=4e3;
K_xy=4e3;
K_yz=2.4e4;

v_mx=0.07;
d_mx=0.005;
b_mx=5e-4;
v_my=0.03;
d_my=0.002;
b_my=6.4e-4;
v_mz=0.072;
d_mz=0.005;
b_mz=4e-4;

v_px=0.128;
d_px=1.6e-4;
v_py=0.128;
d_py=1e-4;
v_pz=0.4;
d_pz=4e-5;

% 物质与时间列表
m_X=[0,init_m_X];
m_Y=[0,init_m_Y];
m_Z=[0,init_m_Z];
X=[0,init_X];
Y=[0,init_Y];
Z=[0,init_Z];

for i=1:4e6
    % 每一步
    a_1=b_mx;
    a_2=v_mx*(Y(end)/K_yx)^n/(1+(Y(end)/K_yx)^n+(Z(end)/K_zx)^n);
    a_3=m_X(end)*d_mx;
    a_4=b_my;
    a_5=v_my*(X(end)/K_xy)^n/(1+(X(end)/K_xy)^n);
    a_6=m_Y(end)*d_my;
    a_7=b_mz;
    a_8=v_mz*(Y(end)/K_yz)^n/(1+(Y(end)/K_yz)^n);
    a_9=m_Z(end)*d_mz;
    a_10=v_px*m_X(end);
    a_11=d_px*X(end);
    a_12=v_py*m_Y(end);
    a_13=d_py*Y(end);
    a_14=v_pz*m_Z(end);
    a_15=d_pz*Z(end);
    a_0=a_1+a_2+a_3+a_4+a_5+a_6+a_7+a_8+a_9+a_10+a_11+a_12+a_13+a_14+a_15;
    
    p_1=rand;
    tau=-log(p_1)/a_0;
    t=t+tau;
    
    s_2=a_1+a_2;
    s_3=s_2+a_3;
    s_5=s_3+a_4+a_5;
    s_6=s_5+a_6;
    s_8=s_6+a_7+a_8;
    s_9=s_8+a_9;
    s_10=s_9+a_10;
    s_11=s_10+a_11;
    s_12=s_11+a_12;
    s_13=s_12+a_13;
    s_14=s_13+a_14;
    s_15=s_14+a_15;
    
    p_2=rand;
    num=p_2*a_0;
    if num<s_2
        m_X=[m_X;[t,m_X(end)+1]];
    elseif num<s_3
        m_X=[m_X;[t,m_X(end)-1]];
    elseif num<s_5
        m_Y=[m_Y;[t,m_Y(end)+1]];
    elseif num<s_6
        m_Y=[m_Y;[t,m_Y(end)-1]];
    elseif num<s_8
        m_Z=[m_Z;[t,m_Z(end)+1]];
    elseif num<s_9
        m_Z=[m_Z;[t,m_Z(end)-1]];
    elseif num<s_10
        X=[X;[t,X(end)+1]];
    elseif num<s_11
        X=[X;[t,X(end)-1]];
    elseif num<s_12
        Y=[Y;[t,Y(end)+1]];
    elseif num<s_13
        Y=[Y;[t,Y(end)-1]];
    elseif num<s_14
        Z=[Z;[t,Z(end)+1]];
    elseif num<s_15
        Z=[Z;[t,Z(end)-1]];
    end
    if t>1e6
        break;
    end
end
%%
plot(X(:,1),X(:,2),'linewidth',2','Color','#0072BD')
hold on
plot(Y(:,1),Y(:,2),'linewidth',2,'Color','#77AC30')
hold on
plot(Z(:,1),Z(:,2),'linewidth',2,'Color','#7E2F8E')
xlabel('Time')
ylabel('Number of Molecules')
xlim([0 1e6])
ylim([0 3e4])
set(gca,'FontSize',16)
legend('[X]','[Y]','[Z]')

toc