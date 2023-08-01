clf

V_y1=readmatrix('aaa.xlsx','Sheet','1p','Range','A1:A1000');
X=readmatrix('aaa.xlsx','Sheet','1p','Range','B1:B1000');
kappa=V_y1.*5;
x=X;

IDsolution_HB=[266,510];
IDsolution_SN=[301,388];

plot(kappa(1:IDsolution_HB(1)),x(1:IDsolution_HB(1)),'k','linewidth',2.5)
hold on
plot(kappa(IDsolution_HB(1):IDsolution_HB(2)),x(IDsolution_HB(1):IDsolution_HB(2)),'--k','linewidth',2.5)
hold on
plot(kappa(IDsolution_HB(2):end),x(IDsolution_HB(2):end),'k','linewidth',2.5)
hold on
plot(kappa(IDsolution_HB(1)),x(IDsolution_HB(1)),'Color','r','Marker','.','MarkerSize',32)
hold on
plot(kappa(IDsolution_SN(1)),x(IDsolution_SN(1)),'Color','b','Marker','.','MarkerSize',24)
hold on
plot(kappa(IDsolution_SN(2)),x(IDsolution_SN(2)),'Color','b','Marker','.','MarkerSize',32)
hold on
plot(kappa(IDsolution_HB(2)),x(IDsolution_HB(2)),'Color','r','Marker','.','MarkerSize',32)
title('1-Parameter Diagram')
set(gca,'FontSize',16)
text(14,0.1,'HB 1','FontSize',12)
text(15,0.2,'SN 1','FontSize',12)
text(12.2,0.45,'SN 2','FontSize',12)
text(17,1.1,'HB 2','FontSize',12)

xlim([10 20])
ylim([0 1.5])
xlabel('\kappa')
ylabel('[x]')