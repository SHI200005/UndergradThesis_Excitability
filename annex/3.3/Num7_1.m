% ���˷��� motif ���� Num5
% �ò����滻�ķ�ʽ���룬���������ʽ�ӱ���
% y=@(x)kappa.*x.^n./(1+x.^n)+beta_y;
% z=@(x)alpha.*(y./K).^n./(1+(y./K).^n)+beta_z;
%
% eta=@(x)(x-beta_x).*(1+y.^n+z.^n)./z.^n;
% b=@(x)-(x-beta_x).*n.*y.^(n-1)./(1+y.^n+z.^n);
% c=@(x)(x-beta_x).*n./z./(1+y.^n+z.^n);
% d=@(x)kappa.*n.*x.^(n-1)./theta_y./(1+x.^n).^2;
% h=@(x)alpha.*n.*(y./K).^(n-1)./theta_z./K./(1+(y./K).^n).^2;
%
% eta=@(x)(x-beta_x).*(1+(kappa.*x.^n./(1+x.^n)+beta_y).^n+(alpha.*((kappa.*x.^n./(1+x.^n)+beta_y)./K).^n./(1+((kappa.*x.^n./(1+x.^n)+beta_y)./K).^n)+beta_z).^n)...,
%     ./(alpha.*((kappa.*x.^n./(1+x.^n)+beta_y)./K).^n./(1+((kappa.*x.^n./(1+x.^n)+beta_y)./K).^n)+beta_z).^n;
% b=@(x)-(x-beta_x).*n.*(kappa.*x.^n./(1+x.^n)+beta_y).^(n-1)./(1+(kappa.*x.^n./(1+x.^n)+beta_y).^n+(alpha.*((kappa.*x.^n./(1+x.^n)+beta_y)./K).^n...,
%     ./(1+((kappa.*x.^n./(1+x.^n)+beta_y)./K).^n)+beta_z).^n);
% c=@(x)(x-beta_x).*n./(alpha.*((kappa.*x.^n./(1+x.^n)+beta_y)./K).^n./(1+((kappa.*x.^n./(1+x.^n)+beta_y)./K).^n)+beta_z)./(1+(kappa.*x.^n./(1+x.^n)+beta_y).^n...,
%     +(alpha.*((kappa.*x.^n./(1+x.^n)+beta_y)./K).^n./(1+((kappa.*x.^n./(1+x.^n)+beta_y)./K).^n)+beta_z).^n);
% d=@(x)kappa.*n.*x.^(n-1)./theta_y./(1+x.^n).^2;
% h=@(x)alpha.*n.*((kappa.*x.^n./(1+x.^n)+beta_y)./K).^(n-1)./theta_z./K./(1+((kappa.*x.^n./(1+x.^n)+beta_y)./K).^n).^2;

tic % ��ʼ��ʱ
paralist=[]; % ���ڴ洢������

% ������ 0.01-100 ��� logspace ȡֵ
for i=1:1e5
    para=Motif('theta_y',1e-2.*10.^(4.*rand),'theta_z',1e-2.*10.^(4.*rand),'K',1e-2.*10.^(4.*rand),'kappa',1e-2.*10.^(4.*rand),...,
        'alpha',1e-2.*10.^(4.*rand),'beta_x',1e-2.*10.^(4.*rand),'beta_y',1e-2.*10.^(4.*rand),'beta_z',1e-2.*10.^(4.*rand));
    paralist=[paralist;para];
end
toc

function [parameter]=Motif(varargin)

% Default Parameter
n=2;% �̶�����
theta_y=1;
theta_z=10;
K=1;
% eta=0.01-100 ������Χ
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
b=@(x)-(x-beta_x).*n.*(kappa.*x.^n./(1+x.^n)+beta_y).^(n-1)./(1+(kappa.*x.^n./(1+x.^n)+beta_y).^n+(alpha.*((kappa.*x.^n./(1+x.^n)+beta_y)./K).^n...,
    ./(1+((kappa.*x.^n./(1+x.^n)+beta_y)./K).^n)+beta_z).^n);
c=@(x)(x-beta_x).*n./(alpha.*((kappa.*x.^n./(1+x.^n)+beta_y)./K).^n./(1+((kappa.*x.^n./(1+x.^n)+beta_y)./K).^n)+beta_z)./(1+(kappa.*x.^n./(1+x.^n)+beta_y).^n...,
    +(alpha.*((kappa.*x.^n./(1+x.^n)+beta_y)./K).^n./(1+((kappa.*x.^n./(1+x.^n)+beta_y)./K).^n)+beta_z).^n);
d=@(x)kappa.*n.*x.^(n-1)./theta_y./(1+x.^n).^2;
h=@(x)alpha.*n.*((kappa.*x.^n./(1+x.^n)+beta_y)./K).^(n-1)./theta_z./K./(1+((kappa.*x.^n./(1+x.^n)+beta_y)./K).^n).^2;

%eigen value ���̵���
A_C=1+1./theta_y+1./theta_z;
B_C=@(x)1./theta_y+1./theta_z+1./theta_y./theta_z-b(x).*d(x);
C_C=@(x)1./theta_y./theta_z-c(x).*d(x).*h(x)-b(x).*d(x)./theta_z;

%�����GSN���Ǳʼ�������ֵ�����е�C_C
GSN=@(x)C_C(x);
%�����GHP��Ӧ�ʼ�������ֵ���̶�ӦA_C*B_C-C_C
GHB=@(x)A_C.*B_C(x)-C_C(x);

% ����etaȡ0-100
% eta.*z.^n./(1+y.^n+z.^n)+beta_x-x;
Naccurate=1e3;%����eta�ľ���
fun_max=@(x)1e2.*(alpha.*((kappa.*x.^n./(1+x.^n)+beta_y)./K).^n./(1+((kappa.*x.^n./(1+x.^n)+beta_y)./K).^n)+beta_z).^n./(1+(kappa.*x.^n./(1+x.^n)+beta_y).^n...,
    +(alpha.*((kappa.*x.^n./(1+x.^n)+beta_y)./K).^n./(1+((kappa.*x.^n./(1+x.^n)+beta_y)./K).^n)+beta_z).^n)+beta_x-x;
x_max=fzero(fun_max,0.1);%��ʼֵһ��ҪȡС
x=linspace(beta_x,x_max,Naccurate);

eta=@(x)(x-beta_x).*(1+(kappa.*x.^n./(1+x.^n)+beta_y).^n+(alpha.*((kappa.*x.^n./(1+x.^n)+beta_y)./K).^n./(1+((kappa.*x.^n./(1+x.^n)+beta_y)./K).^n)+beta_z).^n)...,
    ./(alpha.*((kappa.*x.^n./(1+x.^n)+beta_y)./K).^n./(1+((kappa.*x.^n./(1+x.^n)+beta_y)./K).^n)+beta_z).^n;

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
            if eta(x(IDsolution_HB(1)))>eta(x(IDsolution_SN(2)))
                m=5;
                if IDsolution_HB(1)>=IDsolution_SN(1) %���HB�ֲ���Ǽٵ� Ҳ���Ǳ�Ҫ���������������
                    eta_right=min([eta(x(IDsolution_SN(1))),eta(x(IDsolution_HB(2)))]);
                else % ����HB ���ǵ�һ��
                    eta_right=min([eta(x(IDsolution_HB(1))),eta(x(IDsolution_HB(2)))]);
                end
                eta_left=eta(x(IDsolution_SN(2)));
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
    eta_left=0;
    eta_right=0;
end
if n_SN~=2
    IDsolution_SN=[0,0];% ûɶ�� ֻ�����̽һ�·ֲ������λ�ù�ϵ
end
if n_HB~=2
    IDsolution_HB=[0,0];% ûɶ�� ֻ�����̽һ�·ֲ������λ�ù�ϵ
end

% ����������� ���漸��������ʵûɶ�� ���ǿ��Կ���tri-stability �Լ��ж�Ϊm=0ʱ�ֲ����� �Զ����ԾͿ�
parameter=[m,theta_y,theta_z,K,kappa,alpha,beta_x,beta_y,beta_z,eta_left,eta_right,n_SN,n_HB,IDsolution_HB(1),IDsolution_SN(1),IDsolution_SN(2),IDsolution_HB(2)];
end