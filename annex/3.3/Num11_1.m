% 可兴奋性 motif 搜索 Num12_1
% %用查找替换的方式带入，以免出错，留式子备用
% y=@(x)kappa./(1+x.^n)+beta_y;
% z=@(x)alpha./(1+(x./K).^n)+beta_z;
%
% b=@(x)-(x-beta_x).*n.*y.^(n-1)./(1+y.^n+z.^n);
% c=@(x)(x-beta_x).*n.*(1+y.^n)./z./(1+y.^n+z.^n);
% d=@(x)-kappa.*n.*x.^(n-1)/theta_y./(1+x.^n).^2;
% g=@(x)-alpha.*n.*(x./K).^(n-1)./theta_z./K./(1+(x./K).^n).^2;
% eta=@(x)(x-beta_x).*(1+y.^n+z.^n)./z.*n;
% 
% b=@(x)-(x-beta_x).*n.*(kappa./(1+x.^n)+beta_y).^(n-1)./(1+(kappa./(1+x.^n)+beta_y).^n+(alpha./(1+(x./K).^n)+beta_z).^n);
% c=@(x)(x-beta_x).*n.*(1+(kappa./(1+x.^n)+beta_y).^n)./(alpha./(1+(x./K).^n)+beta_z)./(1+(kappa./(1+x.^n)+beta_y).^n+(alpha./(1+(x./K).^n)+beta_z).^n);
% d=@(x)-kappa.*n.*x.^(n-1)/theta_y./(1+x.^n).^2;
% g=@(x)-alpha.*n.*(x./K).^(n-1)./theta_z./K./(1+(x./K).^n).^2;
% eta=@(x)(x-beta_x).*(1+(kappa./(1+x.^n)+beta_y).^n+(alpha./(1+(x./K).^n)+beta_z).^n)./(alpha./(1+(x./K).^n)+beta_z).^n;

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

%Default Parameter
n=2;%固定不变
theta_y=1;
theta_z=10;
K=1;
%eta=0.01-100;参数范围
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
b=@(x)-(x-beta_x).*n.*(kappa./(1+x.^n)+beta_y).^(n-1)./(1+(kappa./(1+x.^n)+beta_y).^n+(alpha./(1+(x./K).^n)+beta_z).^n);
c=@(x)(x-beta_x).*n.*(1+(kappa./(1+x.^n)+beta_y).^n)./(alpha./(1+(x./K).^n)+beta_z)./(1+(kappa./(1+x.^n)+beta_y).^n+(alpha./(1+(x./K).^n)+beta_z).^n);
d=@(x)-kappa.*n.*x.^(n-1)/theta_y./(1+x.^n).^2;
g=@(x)-alpha.*n.*(x./K).^(n-1)./theta_z./K./(1+(x./K).^n).^2;
%eigen value 方程的项
A_C=1+1./theta_y+1./theta_z;
B_C=@(x)1./theta_y+1./theta_z+1./theta_y./theta_z-c(x).*g(x)-b(x).*d(x);
C_C=@(x)1./theta_y./theta_z-c(x).*g(x)./theta_y-b(x).*d(x)./theta_z;

%这里的GSN就是笔记中特征值方程中的C_C
GSN=@(x)C_C(x);
%这里的GHP对应笔记中特征值方程对应A_C*B_C-C_C
GHB=@(x)A_C.*B_C(x)-C_C(x);

% 参数eta取0-100
% eta.*z.^n./(1+y.^n+z.^n)+beta_x-x
Naccurate=1e3;% 参数eta的精度
fun=@(x)1e2.*(alpha./(1+(x./K).^n)+beta_z).^n./(1+(kappa./(1+x.^n)+beta_y).^n+(alpha./(1+(x./K).^n)+beta_z).^n)+beta_x-x;
x_max=fzero(fun,0.1);
x=linspace(beta_x,x_max,Naccurate);

eta=@(x)(x-beta_x).*(1+(kappa./(1+x.^n)+beta_y).^n+(alpha./(1+(x./K).^n)+beta_z).^n)./(alpha./(1+(x./K).^n)+beta_z).^n;

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
        if IDsolution_HB(2)-IDsolution_SN(2)>1 && IDsolution_HB(1)<IDsolution_SN(2)
            if eta(x(IDsolution_HB(1)))>eta(x(IDsolution_SN(2)))
                m=5;
                if IDsolution_HB(1)>=IDsolution_SN(1) %如果HB分岔点是假的 也就是必要不充分条件的误判
                    eta_right=min([eta(x(IDsolution_SN(1))),eta(x(IDsolution_HB(2)))]);
                else
                    eta_right=min([eta(x(IDsolution_HB(1))),eta(x(IDsolution_HB(2)))]);
                end
                eta_left=eta(x(IDsolution_SN(2)));
            else
                m=4;
            end
        elseif IDsolution_HB(1)<IDsolution_SN(2)
            m=2;
        else
            m=0;
        end
    end
end
% ----此节判定通用！尾

% 只是让输出规整一些 不需要的用0补齐就好
if m~=5
    IDsolution_SN=[0,0];% 没啥用 只是想窥探一下分岔点的相对位置关系
    IDsolution_HB=[0,0];% 没啥用 只是想窥探一下分岔点的相对位置关系
    eta_left=0;
    eta_right=0;
end

% 将参数集输出 后面几个参数其实没啥用 但是可以看看tri-stability 以及判定为m=0时分岔点情况 自动忽略就可
parameter=[m,theta_y,theta_z,K,kappa,alpha,beta_x,beta_y,beta_z,eta_left,eta_right,n_SN,n_HB,IDsolution_HB(1),IDsolution_SN(1),IDsolution_SN(2),IDsolution_HB(2)];
end