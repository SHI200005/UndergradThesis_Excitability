% 可兴奋性 motif 搜索 Num4
% 用查找替换的方式带入，以免出错，留式子备用
% y=@(x)kappa./(1+x.^n)+beta_y;
% z=@(x)alpha./(1+(y./K).^n)+beta_z;
%
% z=@(x)alpha./(1+((kappa./(1+x.^n)+beta_y)./K).^n)+beta_z;
%
% eta=@(x)(x-beta_x).*(1+y.^n+z.^n);
% b=@(x)-(x-beta_x).*n.*y.^(n-1)./(1+y.^n+z.^n);
% c=@(x)-(x-beta_x).*n.*z.^(n-1)./(1+y.^n+z.^n);
% d=@(x)-kappa.*n.*x.^(n-1)./theta_y./(1+x.^n).^2;
% h=@(x)-alpha.*n.*(y./K).^(n-1)./theta_z./K./(1+(y./K).^n).^2;
%
% eta=@(x)(x-beta_x).*(1+(kappa./(1+x.^n)+beta_y).^n+(alpha./(1+((kappa./(1+x.^n)+beta_y)./K).^n)+beta_z).^n);
% b=@(x)-(x-beta_x).*n.*(kappa./(1+x.^n)+beta_y).^(n-1)./(1+(kappa./(1+x.^n)+beta_y).^n+(alpha./(1+((kappa./(1+x.^n)+beta_y)./K).^n)+beta_z).^n);
% c=@(x)-(x-beta_x).*n.*(alpha./(1+((kappa./(1+x.^n)+beta_y)./K).^n)+beta_z).^(n-1)./(1+(kappa./(1+x.^n)+beta_y).^n+(alpha./(1+((kappa./(1+x.^n)+beta_y)./K).^n)+beta_z).^n);
% d=@(x)-kappa.*n.*x.^(n-1)./theta_y./(1+x.^n).^2;
% h=@(x)-alpha.*n.*((kappa./(1+x.^n)+beta_y)./K).^(n-1)./theta_z./K./(1+((kappa./(1+x.^n)+beta_y)./K).^n).^2;

tic % 开始计时
paralist=[]; % 用于存储参数集

% 参数按 0.01-100 随机 logspace 取值
for i=1:1e5
    para=Motif('theta_y',1e-2.*10.^(4.*rand),'theta_z',1e-2.*10.^(4.*rand),'K',1e-2.*10.^(4.*rand),'kappa',1e-2.*10.^(4.*rand),...,
        'alpha',1e-2.*10.^(4.*rand),'beta_x',1e-2.*10.^(4.*rand),'beta_y',1e-2.*10.^(4.*rand),'beta_z',1e-2.*10.^(4.*rand));
    paralist=[paralist;para];
end
toc

function [parameter]=Motif(varargin)

% Default Parameter
n=2;% 固定不变
theta_y=1;
theta_z=10;
K=1;
% eta=0.01-100 参数范围
kappa=14;
alpha=25;
beta_x=0.05;
beta_y=0.5;
beta_z=5;

if ~isempty(varargin)
    for ii=1:length(varargin)/2
        eval([varargin{ii*2-1},'=varargin{ii*2};']);%对参数进行更新
    end
end

%对应笔记中Jacobian矩阵元的符号
b=@(x)-(x-beta_x).*n.*(kappa./(1+x.^n)+beta_y).^(n-1)./(1+(kappa./(1+x.^n)+beta_y).^n+(alpha./(1+((kappa./(1+x.^n)+beta_y)./K).^n)+beta_z).^n);
c=@(x)-(x-beta_x).*n.*(alpha./(1+((kappa./(1+x.^n)+beta_y)./K).^n)+beta_z).^(n-1)./(1+(kappa./(1+x.^n)+beta_y).^n+(alpha./(1+((kappa./(1+x.^n)+beta_y)./K).^n)+beta_z).^n);
d=@(x)-kappa.*n.*x.^(n-1)./theta_y./(1+x.^n).^2;
h=@(x)-alpha.*n.*((kappa./(1+x.^n)+beta_y)./K).^(n-1)./theta_z./K./(1+((kappa./(1+x.^n)+beta_y)./K).^n).^2;

%eigen value 方程的项
A_C=1+1./theta_y+1./theta_z;
B_C=@(x)1./theta_y+1./theta_z+1./theta_y./theta_z-b(x).*d(x);
C_C=@(x)1./theta_y./theta_z-c(x).*d(x).*h(x)-b(x).*d(x)./theta_z;

%这里的GSN就是笔记中特征值方程中的C_C
GSN=@(x)C_C(x);
%这里的GHP对应笔记中特征值方程对应A_C*B_C-C_C
GHB=@(x)A_C.*B_C(x)-C_C(x);

% 参数eta取0-100
% eta./(1+y.^n+z.^n)+beta_x-x;
Naccurate=1e3;%参数eta的精度
fun_max=@(x)1e2./(1+(kappa./(1+x.^n)+beta_y).^n+(alpha./(1+((kappa./(1+x.^n)+beta_y)./K).^n)+beta_z).^n)+beta_x-x;
x_max=fzero(fun_max,0.1);%初始值一定要取小
x=linspace(beta_x,x_max,Naccurate);

eta=@(x)(x-beta_x).*(1+(kappa./(1+x.^n)+beta_y).^n+(alpha./(1+((kappa./(1+x.^n)+beta_y)./K).^n)+beta_z).^n);

% ----此节判定通用！头
%先扫描SN分岔
temp1=GSN(x);
IDsolution_SN=find((temp1(1:end-1).*temp1(2:end))<0);
n_SN=size(IDsolution_SN,2);
%再扫描HB分岔
temp2=GHB(x);
IDsolution_HB=find((temp2(1:end-1).*temp2(2:end))<0);
n_HB=size(IDsolution_HB,2);

m=0; % m 代表行为的类别 如果所有情况都不符合，则暂不考虑
if n_SN==0
    if n_HB==0
        m=1;%单稳1
    elseif n_HB==2
        m=3;%振荡3
    end
elseif n_SN==2
    if n_HB==0
        m=2;%双稳
    elseif n_HB==2
        if IDsolution_HB(2)-IDsolution_SN(2)>1 && IDsolution_HB(1)<IDsolution_SN(2) % 别是俩HB都在外边
            if eta(x(IDsolution_HB(1)))>eta(x(IDsolution_SN(2)))
                m=5;
                if IDsolution_HB(1)>=IDsolution_SN(1) %如果HB分岔点是假的 也就是必要不充分条件的误判
                    eta_right=min([eta(x(IDsolution_SN(1))),eta(x(IDsolution_HB(2)))]);
                else % 两个HB 考虑第一个
                    eta_right=min([eta(x(IDsolution_HB(1))),eta(x(IDsolution_HB(2)))]);
                end
                eta_left=eta(x(IDsolution_SN(2)));
            else
                m=4; 
            end
        elseif IDsolution_HB(1)<IDsolution_SN(2) % 并不是两个都在外面 但是可能HB2和SN2在一起 因为变号点可能是左边界也可能是右边界 一起变号不一定分开了
            m=2; % HB2和SN2靠在一起就是双稳
        else
            m=0; % 俩HB都在外面咱就不考虑了喵喵
        end
    end
end
% ----此节判定通用！尾

% 只是让输出规整一些 不需要的用0补齐就好
if m~=5
    eta_left=0;
    eta_right=0;
end
if n_SN~=2
    IDsolution_SN=[0,0];% 没啥用 只是想窥探一下分岔点的相对位置关系
end
if n_HB~=2
    IDsolution_HB=[0,0];% 没啥用 只是想窥探一下分岔点的相对位置关系
end

% 将参数集输出 后面几个参数其实没啥用 但是可以看看tri-stability 以及判定为m=0时分岔点情况 自动忽略就可
parameter=[m,theta_y,theta_z,K,kappa,alpha,beta_x,beta_y,beta_z,eta_left,eta_right,n_SN,n_HB,IDsolution_HB(1),IDsolution_SN(1),IDsolution_SN(2),IDsolution_HB(2)];
end