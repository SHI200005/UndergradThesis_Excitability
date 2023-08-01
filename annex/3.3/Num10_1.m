% ���˷��� motif ���� Num10
% �ò����滻�ķ�ʽ���룬���������ʽ�ӱ���
% y=@(x)kappa.*x.^n./(1+x.^n)+beta_y;
% z=@(x)alpha./(1+(x./K).^n)+beta_z;
%
% b=@(x)n.*y.^(n-1).*(miu.*(1+z.^n)-niu.*z.^n)./(1+y.^n+z.^n).^2 % ��û������ ���
% c=@(x)n.*z.^(n-1).*(-miu.*y.^n+niu.*(1+y.^n))./(1+y.^n+z.^n).^2
% miu=@(x)((1+y.^n+z.^n).*(x-beta_x)-niu.*z.^n)./y.^n
%
% b=@(x)n.*y.^(n-1).*((((1+y.^n+z.^n).*(x-beta_x)-niu.*z.^n)./y.^n).*(1+z.^n)-niu.*z.^n)./(1+y.^n+z.^n).^2
% c=@(x)n.*z.^(n-1).*(-(((1+y.^n+z.^n).*(x-beta_x)-niu.*z.^n)./y.^n).*y.^n+niu.*(1+y.^n))./(1+y.^n+z.^n).^2
% d=@(x)kappa.*n.*x.^(n-1)/theta_y./(1+x.^n).^2;
% g=@(x)-alpha.*n.*(x./K).^(n-1)./theta_z./K./(1+(x./K).^n).^2;
%
% miu=@(x)((1+y.^n+z.^n).*(x-beta_x)-niu.*z.^n)./y.^n
%
% b=@(x)n.*(kappa.*x.^n./(1+x.^n)+beta_y).^(n-1).*((((1+(kappa.*x.^n./(1+x.^n)+beta_y).^n+(alpha./(1+(x./K).^n)+beta_z).^n)...,
%     .*(x-beta_x)-niu.*(alpha./(1+(x./K).^n)+beta_z).^n)./(kappa.*x.^n./(1+x.^n)+beta_y).^n).*(1+(alpha./(1+(x./K).^n)+beta_z).^n)...,
%     -niu.*(alpha./(1+(x./K).^n)+beta_z).^n)./(1+(kappa.*x.^n./(1+x.^n)+beta_y).^n+(alpha./(1+(x./K).^n)+beta_z).^n).^2;
% c=@(x)n.*(alpha./(1+(x./K).^n)+beta_z).^(n-1).*(-(((1+(kappa.*x.^n./(1+x.^n)+beta_y).^n+(alpha./(1+(x./K).^n)+beta_z).^n).*(x-beta_x)...,
%     -niu.*(alpha./(1+(x./K).^n)+beta_z).^n)./(kappa.*x.^n./(1+x.^n)+beta_y).^n).*(kappa.*x.^n./(1+x.^n)+beta_y).^n+niu...,
%     .*(1+(kappa.*x.^n./(1+x.^n)+beta_y).^n))./(1+(kappa.*x.^n./(1+x.^n)+beta_y).^n+(alpha./(1+(x./K).^n)+beta_z).^n).^2;
% d=@(x)kappa.*n.*x.^(n-1)/theta_y./(1+x.^n).^2;
% g=@(x)-alpha.*n.*(x./K).^(n-1)./theta_z./K./(1+(x./K).^n).^2;
%
% miu=@(x)((1+(kappa.*x.^n./(1+x.^n)+beta_y).^n+(alpha./(1+(x./K).^n)+beta_z).^n).*(x-beta_x)-niu.*(alpha./(1+(x./K).^n)+beta_z).^n)./(kappa.*x.^n./(1+x.^n)+beta_y).^n

tic % ��ʼ��ʱ
paralist=[]; % ���ڴ洢������

% ������ 0.01-100 ��� logspace ȡֵ
for i=1:1e5
    para=Motif('niu',1e-2.*10.^(4.*rand),'theta_y',1e-2.*10.^(4.*rand),'theta_z',1e-2.*10.^(4.*rand),'K',1e-2.*10.^(4.*rand),'kappa',1e-2.*10.^(4.*rand),...,
        'alpha',1e-2.*10.^(4.*rand),'beta_x',1e-2.*10.^(4.*rand),'beta_y',1e-2.*10.^(4.*rand),'beta_z',1e-2.*10.^(4.*rand));
    paralist=[paralist;para];
end
toc

function [parameter]=Motif(varargin)

% Default Parameter ��������������Ƿ��Ͽ��˷��Ե�һ������
n=2;% �̶�����
niu=1;
theta_y=1;
theta_z=10;
K=1;
% miu=0.01-100 ϣ���Ĳ�����Χ Ȼ����ʱ����ԸΥ ι ��������ԸΥ�ˣ�
kappa=14;
alpha=25;
beta_x=0.05;
beta_y=0.5;
beta_z=5;

if ~isempty(varargin)
    for ii=1:length(varargin)/2
        eval([varargin{ii*2-1},'=varargin{ii*2};']);%�Բ������и���
    end
end

%��Ӧ�ʼ���Jacobian����Ԫ�ķ���
b=@(x)n.*(kappa.*x.^n./(1+x.^n)+beta_y).^(n-1).*((((1+(kappa.*x.^n./(1+x.^n)+beta_y).^n+(alpha./(1+(x./K).^n)+beta_z).^n)...,
    .*(x-beta_x)-niu.*(alpha./(1+(x./K).^n)+beta_z).^n)./(kappa.*x.^n./(1+x.^n)+beta_y).^n).*(1+(alpha./(1+(x./K).^n)+beta_z).^n)...,
    -niu.*(alpha./(1+(x./K).^n)+beta_z).^n)./(1+(kappa.*x.^n./(1+x.^n)+beta_y).^n+(alpha./(1+(x./K).^n)+beta_z).^n).^2;
c=@(x)n.*(alpha./(1+(x./K).^n)+beta_z).^(n-1).*(-(((1+(kappa.*x.^n./(1+x.^n)+beta_y).^n+(alpha./(1+(x./K).^n)+beta_z).^n).*(x-beta_x)...,
    -niu.*(alpha./(1+(x./K).^n)+beta_z).^n)./(kappa.*x.^n./(1+x.^n)+beta_y).^n).*(kappa.*x.^n./(1+x.^n)+beta_y).^n+niu...,
    .*(1+(kappa.*x.^n./(1+x.^n)+beta_y).^n))./(1+(kappa.*x.^n./(1+x.^n)+beta_y).^n+(alpha./(1+(x./K).^n)+beta_z).^n).^2;
d=@(x)kappa.*n.*x.^(n-1)/theta_y./(1+x.^n).^2;
g=@(x)-alpha.*n.*(x./K).^(n-1)./theta_z./K./(1+(x./K).^n).^2;

%eigen value ���̵���
A_C=1+1./theta_y+1./theta_z;
B_C=@(x)1./theta_y+1./theta_z+1./theta_y./theta_z-c(x).*g(x)-b(x).*d(x);
C_C=@(x)1./theta_y./theta_z-c(x).*g(x)./theta_y-b(x).*d(x)./theta_z;

%�����GSN���Ǳʼ�������ֵ�����е�C_C
GSN=@(x)C_C(x);
%�����GHP��Ӧ�ʼ�������ֵ���̶�ӦA_C*B_C-C_C
GHB=@(x)A_C.*B_C(x)-C_C(x);

% ����miuȡ0-100
% (miu.*y.^n+niu.*z.^n)/(1+y.^n+z.^n)+beta_x-x;
% niu.*z.^n./(1+y.^n+z.^n)+beta_x-x;
Naccurate=1e3;%����eta�ľ���
fun_max=@(x)(1e2.*(kappa.*x.^n./(1+x.^n)+beta_y).^n+niu.*(alpha./(1+(x./K).^n)+beta_z).^n)/(1+(kappa.*x.^n./(1+x.^n)+beta_y).^n+(alpha./(1+(x./K).^n)+beta_z).^n)+beta_x-x;
x_max=fzero(fun_max,0.1);%��ʼֵһ��ҪȡС
fun_min=@(x)niu.*(alpha./(1+(x./K).^n)+beta_z).^n./(1+(kappa.*x.^n./(1+x.^n)+beta_y).^n+(alpha./(1+(x./K).^n)+beta_z).^n)+beta_x-x; % ע��Ŷ ������Сֵ�ɲ���beta_x�� ����
x_min=fzero(fun_min,0.01);
x=linspace(x_min,x_max,Naccurate);

miu=@(x)((1+(kappa.*x.^n./(1+x.^n)+beta_y).^n+(alpha./(1+(x./K).^n)+beta_z).^n).*(x-beta_x)-niu.*(alpha./(1+(x./K).^n)+beta_z).^n)./(kappa.*x.^n./(1+x.^n)+beta_y).^n;

% ----�˽��ж�ͨ�ã�ͷ
%��ɨ��SN�ֲ�
temp1=GSN(x);
IDsolution_SN=find((temp1(1:end-1).*temp1(2:end))<0);
n_SN=size(IDsolution_SN,2);
%��ɨ��HB�ֲ�
temp2=GHB(x);
IDsolution_HB=find((temp2(1:end-1).*temp2(2:end))<0);
n_HB=size(IDsolution_HB,2);

m=0; % m ������Ϊ����� �����������������ϣ����ݲ�����
if n_SN==0
    if n_HB==0
        m=1;%����1
    elseif n_HB==2
        m=3;%��3
    end
elseif n_SN==2
    if n_HB==0
        m=2;%˫��
    elseif n_HB==2
        if IDsolution_HB(2)-IDsolution_SN(2)>1 && IDsolution_HB(1)<IDsolution_SN(2) % ������HB�������
            if miu(x(IDsolution_HB(1)))>miu(x(IDsolution_SN(2)))
                m=5;
                if IDsolution_HB(1)>=IDsolution_SN(1) %���HB�ֲ���Ǽٵ� Ҳ���Ǳ�Ҫ���������������
                    miu_right=min([miu(x(IDsolution_SN(1))),miu(x(IDsolution_HB(2)))]);
                else % ����HB ���ǵ�һ��
                    miu_right=min([miu(x(IDsolution_HB(1))),miu(x(IDsolution_HB(2)))]);
                end
                miu_left=miu(x(IDsolution_SN(2)));
            else
                m=4; 
            end
        elseif IDsolution_HB(1)<IDsolution_SN(2) % ������������������ ���ǿ���HB2��SN2��һ�� ��Ϊ��ŵ��������߽�Ҳ�������ұ߽� һ���Ų�һ���ֿ���
            m=2; % HB2��SN2����һ�����˫��
        else
            m=0; % ��HB���������۾Ͳ�����������
        end
    end
end
% ----�˽��ж�ͨ�ã�β

% ֻ�����������һЩ ����Ҫ����0����ͺ�
if m~=5
    miu_left=0;
    miu_right=0;
end
if n_SN~=2
    IDsolution_SN=[0,0];% ûɶ�� ֻ�����̽һ�·ֲ������λ�ù�ϵ
end
if n_HB~=2
    IDsolution_HB=[0,0];% ûɶ�� ֻ�����̽һ�·ֲ������λ�ù�ϵ
end

% ����������� ���漸��������ʵûɶ�� ���ǿ��Կ���tri-stability �Լ��ж�Ϊm=0ʱ�ֲ����� �Զ����ԾͿ�
parameter=[m,niu,theta_y,theta_z,K,kappa,alpha,beta_x,beta_y,beta_z,miu_left,miu_right,n_SN,n_HB,IDsolution_HB(1),IDsolution_SN(1),IDsolution_SN(2),IDsolution_HB(2)];
end