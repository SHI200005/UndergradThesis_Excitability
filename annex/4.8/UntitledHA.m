clf
HB3=readmatrix('aaa.xlsx','Sheet','HB3','Range','A1:Z2000');
SN=readmatrix('aaa.xlsx','Sheet','SN','Range','A1:Z2000');
HB11=readmatrix('aaa.xlsx','Sheet','HB11','Range','A1:Z2000');

plot(HB3(:,1).*5,HB3(:,2).*50,'Color','#0072BD','linewidth',2.5)
hold on
plot(SN(:,1).*5,SN(:,2).*50,'Color','#D95319','linewidth',2.5)
hold on
plot(HB11(:,1).*5,HB11(:,2).*50,'Color','#7E2F8E','linewidth',2.5)
xlim([0.8.*5 3.8.*5])
ylim([0.*50 1.*50])
set(gca,'FontSize',16)
xlabel('\kappa','FontSize',20)
ylabel('\alpha','FontSize',20)
title('2-Parameter Diagram')
text(15.5,45,'HB 1','FontSize',16)
text(12,5,'SN 1','FontSize',16)
text(11,24,'SN 2','FontSize',16)
text(17,24,'HB 2','FontSize',16)
