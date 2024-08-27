
clear all

% Add package to read NPY files into matlab
addpath('C:/Users/kevin/OneDrive/Documents/NGG_PhD/Matlab_package/npy-matlab')

% Load in prediction results
folder_name='../../Step_2/2_Server_all/BASELINE/baseline_results_031324/all_network';
file_name_pred1_samp1='General_PB1_prediction_testA.npy';
file_name_pred1_samp2='General_PB1_prediction_testB.npy';
pred_PB1_samp1 = readNPY([folder_name '/' file_name_pred1_samp1]);
pred_PB1_samp2 = readNPY([folder_name '/' file_name_pred1_samp2]);


% put predictions on same scale as actual
pred_PB1_samp1 = zscore(pred_PB1_samp1);
pred_PB1_samp2 = zscore(pred_PB1_samp2);

% fourth column has the matched_group
data_for_ridge=xlsread('../../Step_2/2_Server_all/BASELINE/data_for_ridge_pfac_BASELINE.csv')

% Load in actual pfac scores
actual_PB1 = data_for_ridge(:,3);

% split data out by matched group
actual_PB1_samp1 = actual_PB1(data_for_ridge(:,4)==1);
actual_PB1_samp2 = actual_PB1(data_for_ridge(:,4)==2);


%% Plotting

% load yellow colors
color_fill_1 = [0.78, 0.60, 0.01]
color_fill_2 = [1, 0.75, 0]

% PB1
figure()
scatter(actual_PB1_samp2,pred_PB1_samp2,15,color_fill_2,'Marker','^','MarkerFaceColor',color_fill_2)
[r,p,ci,stats]=corrcoef(actual_PB1_samp2,pred_PB1_samp2)
PB1_samp1_r = r(1,2)
PB1_samp1_p = p(1,2)
PB1_samp1_CI = [ci(1,2),stats(1,2)]
hold on
scatter(actual_PB1_samp1,pred_PB1_samp1,75,color_fill_1,'Marker','.')
[r,p,ci,stats]=corrcoef(actual_PB1_samp1,pred_PB1_samp1);
PB1_samp2_r = r(1,2)
PB1_samp2_p = p(1,2)
PB1_samp2_CI = [ci(1,2),stats(1,2)]
hold on
linearFit = polyfit(actual_PB1_samp1,pred_PB1_samp1,1);
hline = refline(linearFit);
hline.Color=color_fill_1;
hline.LineWidth=2;
hold on
linearFit = polyfit(actual_PB1_samp2,pred_PB1_samp2,1);
hline = refline(linearFit);
hline.Color=color_fill_2;
hline.LineWidth=2;
title("P-factor Model Performance",'FontSize',20)
xlabel("Actual P-factor Score",'FontSize',16)
ylabel("Model-Fit P-factor Score",'FontSize',16)
axis([-3 4 -4 4])

