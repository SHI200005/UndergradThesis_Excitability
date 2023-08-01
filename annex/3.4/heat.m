clf

heat_1=readmatrix('1.xlsx','Sheet','1','Range','J1:Q8');
heat_3=readmatrix('1.xlsx','Sheet','3','Range','J1:Q8');
heat_4=readmatrix('1.xlsx','Sheet','4','Range','J1:Q8');
heat_9=readmatrix('1.xlsx','Sheet','9','Range','J1:Q8');
heat_11=readmatrix('1.xlsx','Sheet','11','Range','J1:Q8');
heat_12=readmatrix('1.xlsx','Sheet','12','Range','J1:Q8');

xname = {'\theta_y','\theta_z','K','\kappa','\alpha','\beta_x','\beta_y','\beta_z'};
yname = {'31.60~100.0','10.00~31.60','3.160~10.00','1.000~3.160','0.316~1.000','0.100~0.316','0.0316~0.100','0.010~0.0316'};
%%
subplot(1,2,1)
h_1=heatmap(xname,yname,heat_1);
h_1.CellLabelColor = 'none';
h_1.Title='#1';
% h_1.FontSize=14;

subplot(1,2,2)
h_9=heatmap(xname,yname,heat_9);
h_9.CellLabelColor = 'none';
h_9.Title='#9';
% h_9.FontSize=14;
%%
% subplot(1,2,1)
% h_3=heatmap(xname,yname,heat_3);
% h_3.CellLabelColor = 'none';
% h_3.Title='#3';
% 
% subplot(1,2,2)
% h_11=heatmap(xname,yname,heat_11);
% h_11.CellLabelColor = 'none';
% h_11.Title='#11';

% subplot(1,2,1)
% h_4=heatmap(xname,yname,heat_4);
% h_4.CellLabelColor = 'none';
% h_4.Title='#4';
% 
% subplot(1,2,2)
% h_12=heatmap(xname,yname,heat_12);
% h_12.CellLabelColor = 'none';
% h_12.Title='#12';