clear
clf
tic
outputlist=[];
robustlist=[];
list=readmatrix('11.xlsx','Sheet','1','Range','A1:Z200');

for i=1:size(list,1)
    
    outputlist=[];
    theta_y=list(i,1);
    theta_z=list(i,2);
    K=list(i,3);
    kappa=list(i,4);
    alpha=list(i,5);
    beta_x=list(i,6);
    beta_y=list(i,7);
    beta_z=list(i,8);
    x_max=list(i,9);
    for ii=1:8 % 每组数据8个参数
        
        input=[theta_y,theta_z,K,kappa,alpha,beta_x,beta_y,beta_z,x_max,ii];
        output_min_max=check(input);
        outputlist=[outputlist,output_min_max]; % 横着来16个
    end
    robustlist=[robustlist;outputlist]; % 竖着来几十个
end
%%
x = 1:8*size(list,1);
for i=1:8
    y_min=robustlist(:,2*i-1);
    y_max=robustlist(:,2*i);
    y=[y_min,y_max-y_min];
    b=bar(x(1+size(list,1)*(i-1):size(list,1)*i),y,'stacked');
    set(gca,'yscal','log')
    set(b,'edgecolor','none')
    set(b(1),'facecolor','white')
    hold on
end
title('#12','FontSize',16)
set(gca,'FontSize',16,'xtick',[],'xticklabel',[])

%%
toc

function [output_min_max]=check(input)
n=2;
ii=input(10);
temp=input(ii);
for iii=linspace(-1,-40,40) %
    output=0;
    input(ii)=temp.*10.^(iii.*0.1);
    theta_y=input(1);
    theta_z=input(2);
    K=input(3);
    kappa=input(4);
    alpha=input(5);
    beta_x=input(6);
    beta_y=input(7);
    beta_z=input(8);
    x_max=input(9);
    x=linspace(beta_x,x_max,1e3);
    % 对应笔记中Jacobian矩阵元的符号
    % ---------------------------------------------
    b=@(x)-(x-beta_x).*n.*(kappa./(1+x.^n)+beta_y).^(n-1)./(1+(kappa./(1+x.^n)+beta_y).^n+(alpha.*(x./K).^n./(1+(x./K).^n)+beta_z).^n);
    c=@(x)-(x-beta_x).*n.*(alpha.*(x./K).^n./(1+(x./K).^n)+beta_z).^(n-1)./(1+(kappa./(1+x.^n)+beta_y).^n+(alpha.*(x./K).^n./(1+(x./K).^n)+beta_z).^n);
    d=@(x)-kappa.*n.*x.^(n-1)/theta_y./(1+x.^n).^2;
    g=@(x)alpha.*n.*(x./K).^(n-1)./theta_z./K./(1+(x./K).^n).^2;
    % 特征值方程中
    A_C=1+1./theta_y+1./theta_z;
    B_C=@(x)1./theta_y+1./theta_z+1./theta_y./theta_z-c(x).*g(x)-b(x).*d(x);
    C_C=@(x)1./theta_y./theta_z-c(x).*g(x)./theta_y-b(x).*d(x)./theta_z;
    % ---------------------------------------------
    %这里的GSN就是笔记中特征值方程中的C_C
    GSN=@(x)C_C(x);
    %这里的GHP对应笔记中特征值方程对应A_C*B_C-C_C
    GHB=@(x)A_C.*B_C(x)-C_C(x);
    temp1=GSN(x);
    IDsolution_SN=find((temp1(1:end-1).*temp1(2:end))<0);
    n_SN=size(IDsolution_SN,2);
    %再扫描HB分岔
    temp2=GHB(x);
    IDsolution_HB=find((temp2(1:end-1).*temp2(2:end))<0);
    n_HB=size(IDsolution_HB,2);
    if n_SN~=2 || n_HB~=2
        output=input(ii).*10.^(0.1); %
    else
        if IDsolution_HB(2)<IDsolution_SN(2)+1
            output=input(ii).*10.^(0.1);
        else
            eta=@(x)(x-beta_x).*(1+(kappa./(1+x.^n)+beta_y).^n+(alpha.*(x./K).^n./(1+(x./K).^n)+beta_z).^n);
            if eta(x(IDsolution_HB(1)))<eta(x(IDsolution_SN(2)))
                output=input(ii).*10.^(0.1);
            end
        end
    end
    if output~=0
        break;
    end
    if input(ii)<0.01
        output=0.01;
        break;
    end
end
output_min=output;

for iii=linspace(1,40,40) %
    output=0;
    input(ii)=temp.*10.^(iii.*0.1);
    theta_y=input(1);
    theta_z=input(2);
    K=input(3);
    kappa=input(4);
    alpha=input(5);
    beta_x=input(6);
    beta_y=input(7);
    beta_z=input(8);
    x_max=input(9);
    x=linspace(beta_x,x_max,1e3);
    % 对应笔记中Jacobian矩阵元的符号
    b=@(x)-(x-beta_x).*n.*(kappa./(1+x.^n)+beta_y).^(n-1)./(1+(kappa./(1+x.^n)+beta_y).^n+(alpha.*(x./K).^n./(1+(x./K).^n)+beta_z).^n);
    c=@(x)-(x-beta_x).*n.*(alpha.*(x./K).^n./(1+(x./K).^n)+beta_z).^(n-1)./(1+(kappa./(1+x.^n)+beta_y).^n+(alpha.*(x./K).^n./(1+(x./K).^n)+beta_z).^n);
    d=@(x)-kappa.*n.*x.^(n-1)/theta_y./(1+x.^n).^2;
    g=@(x)alpha.*n.*(x./K).^(n-1)./theta_z./K./(1+(x./K).^n).^2;
    % 特征值方程中
    A_C=1+1./theta_y+1./theta_z;
    B_C=@(x)1./theta_y+1./theta_z+1./theta_y./theta_z-c(x).*g(x)-b(x).*d(x);
    C_C=@(x)1./theta_y./theta_z-c(x).*g(x)./theta_y-b(x).*d(x)./theta_z;
    %这里的GSN就是笔记中特征值方程中的C_C
    GSN=@(x)C_C(x);
    %这里的GHP对应笔记中特征值方程对应A_C*B_C-C_C
    GHB=@(x)A_C.*B_C(x)-C_C(x);
    temp1=GSN(x);
    IDsolution_SN=find((temp1(1:end-1).*temp1(2:end))<0);
    n_SN=size(IDsolution_SN,2);
    %再扫描HB分岔
    temp2=GHB(x);
    IDsolution_HB=find((temp2(1:end-1).*temp2(2:end))<0);
    n_HB=size(IDsolution_HB,2);
    if n_SN~=2 || n_HB~=2
        output=input(ii).*10.^(-0.1); %
    else
        if IDsolution_HB(2)<IDsolution_SN(2)+1
            output=input(ii).*10.^(-0.1);
        else
            eta=@(x)(x-beta_x).*(1+(kappa./(1+x.^n)+beta_y).^n+(alpha.*(x./K).^n./(1+(x./K).^n)+beta_z).^n);
            if eta(x(IDsolution_HB(1)))<eta(x(IDsolution_SN(2)))
                output=input(ii).*10.^(-0.1);
            end
        end
    end
    if output~=0
        break;
    end
    if input(ii)>100
        output=100;
        break;
    end
end
output_min_max=[output_min,output];
end