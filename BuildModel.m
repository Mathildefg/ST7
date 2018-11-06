clc; clear; close all;
%% ST7 - Gaussian Process Regression (GPR)
%This script uses data collected using the MATLAB-app
%'CollectTrainingData', to build a GPR-model.

%% Loading and structuring data
% Prompt the user to choose the folder containing the data of the
% respective test subject. 
fileFolder = dir(uigetdir);
fileFolder = fileFolder.folder;

% Get a list of all .mat-files in the folder.
listOfFiles = dir(fullfile(fileFolder, '*.mat'));
% Load the files 
EMG= []; RMS= []; GeneratedProfile = []; movement=[];
for k = 1 : length(listOfFiles)

  load(listOfFiles(k).name, 'emg');
  data{k} = emg.values;
  EMG = [EMG; data{1,k}];
  
  RMS_temp = [];
  for i=1 : 8
      RMS_temp(:,i) = transpose(my_rms(data{1,k}(:,i),40,20,0));
     
  end
  RMS = [RMS; RMS_temp];
  
  GeneratedProfile = [GeneratedProfile my_rms(data{1,k}, 40, 20, 0)];
  
  %To be used for indexing - the lengths of the RMS from each movement:
  f(k) = length(RMS_temp);
  
end
clear emg RMS_temp;

GeneratedProfile= transpose(GeneratedProfile);

%Assigning a number for each movement to be used in indexing. Feasible
%because the files are always loaded alphabetically9
movement(1:sum(f(1:3))) =2;                   %Extension = 2
movement(sum(f(1:3))+1: sum(f(1:6)))=1;       %Flexion= 1
movement(sum(f(1:6))+1:sum(f(1:9)))=6;        %Pronation = 6
movement(sum(f(1:9))+1:sum(f(1:12)))=3;       %Radial deviation = 3
movement(sum(f(1:12))+1:sum(f(1:15)))=7;      %Rest = 7
movement(sum(f(1:15))+1:sum(f(1:19)))= 5;     %Supination = 5
movement(sum(f(1:19))+1:sum(f(1:21)))=4;      %Ulnar deviation = 4

%The index vector
[dofs, i, ~] = unique(movement);
i = sort(i);


figure(1)
subplot(3,1,1)
plot(1:length(EMG),EMG)
subplot(3,1,2)
plot(1:length(RMS),RMS)
subplot(3,1,3)
plot(1:length(GeneratedProfile), GeneratedProfile)


%% Arrangement in predictors and target

%Predictor - RMS of each EMG-channel of the myo armband for 25%, 50%, 75% of MVC. 
x_ext = RMS(i(1):i(2)-1,:);  %Vi tr�kker 1 fra fordi i(2) er f�rste sample i n�ste bev�gelse.
x_flex = RMS(i(2):i(3)-1,:);
x_rd = RMS(i(3):i(4)-1,:);
x_ud = RMS(i(4):i(5)-1,:);
x_sup = RMS(i(5):i(6)-1,:);
x_pro = RMS(i(6):i(7)-1,:);
x_rest = RMS(i(7-1):end,:);

% Setting rest to zero
GeneratedProfile(i(7):end,:)= 0;

%Target values - The generated profile. 
y_flex = GeneratedProfile(i(1):i(2)-1,:);  
y_ext = GeneratedProfile(i(2):i(3)-1,:);
y_rd = GeneratedProfile(i(3):i(4)-1,:);
y_ud = GeneratedProfile(i(4):i(5)-1,:);
y_sup = GeneratedProfile(i(5):i(6)-1,:);
y_pro = GeneratedProfile(i(6):i(7)-1,:);
y_rest = GeneratedProfile(i(7):end,:);


%
 y1 = [y_flex;y_rest]; x1 = [y_flex;y_rest];
 y2 = [y_ext;y_rest];  x2 = [y_ext;y_rest];
 y3 = [y_rd;y_rest];   x3 = [y_rd;y_rest];
 y4 = [y_ud;y_rest];   x4 = [y_ud;y_rest];
 y5 = [y_sup;y_rest];  x5 = [y_sup;y_rest];
 y6 = [y_pro;y_rest];  x6 = [y_pro;y_rest];

 


 
%% Gaussian Process Regression - Stort set kopieret fra Baumeister.
optimization = input('Optimization? Y/N [N]: ','s');
if ~strcmp(optimization,'Y') && ~strcmp(optimization,'y')
    optimization = 'N';
elseif strcmp(optimization,'y')
    optimization = 'Y';
end

switch optimization %the little o after models stands for 'old'.
    case 'Y'        % if optimization method is chosen
        % Gaussian Regression - Optimization
        gprMdl_dof1o = fitrgp(x1,y1,'FitMethod','exact','PredictMethod','exact','Basisfunction','none','OptimizeHyperparameters',{'Sigma'},'HyperparameterOptimizationOptions', struct('MaxObjectiveEvaluations',20));
        gprMdl_dof2o = fitrgp(x2,y2,'FitMethod','exact','PredictMethod','exact','Basisfunction','none','OptimizeHyperparameters',{'Sigma'},'HyperparameterOptimizationOptions', struct('MaxObjectiveEvaluations',20));
        gprMdl_dof3o = fitrgp(x3,y3,'FitMethod','exact','PredictMethod','exact','Basisfunction','none','OptimizeHyperparameters',{'Sigma'},'HyperparameterOptimizationOptions', struct('MaxObjectiveEvaluations',20));
        gprMdl_dof4o = fitrgp(x4,y4,'FitMethod','exact','PredictMethod','exact','Basisfunction','none','OptimizeHyperparameters',{'Sigma'},'HyperparameterOptimizationOptions', struct('MaxObjectiveEvaluations',20));
        gprMdl_dof4o = fitrgp(x5,y5,'FitMethod','exact','PredictMethod','exact','Basisfunction','none','OptimizeHyperparameters',{'Sigma'},'HyperparameterOptimizationOptions', struct('MaxObjectiveEvaluations',20));
        gprMdl_dof4o = fitrgp(x6,y6,'FitMethod','exact','PredictMethod','exact','Basisfunction','none','OptimizeHyperparameters',{'Sigma'},'HyperparameterOptimizationOptions', struct('MaxObjectiveEvaluations',20));
        
    case 'N'        % if non-optimization (already known hyperparameters) method is chosen
        % Gaussian Regression - Already known hyperparameters
        gprMdl_dof1o = fitrgp(x1,y1,'Fitmethod','exact','PredictMethod','exact','Sigma',0.06,'KernelParameters',[0.15;0.35],'Basisfunction','none');
        gprMdl_dof2o = fitrgp(x2,y2,'Fitmethod','exact','PredictMethod','exact','Sigma',0.06,'KernelParameters',[0.15;0.35],'Basisfunction','none');
        gprMdl_dof3o = fitrgp(x3,y3,'Fitmethod','exact','PredictMethod','exact','Sigma',0.06,'KernelParameters',[0.15;0.35],'Basisfunction','none');
        gprMdl_dof4o = fitrgp(x4,y4,'Fitmethod','exact','PredictMethod','exact','Sigma',0.06,'KernelParameters',[0.15;0.35],'Basisfunction','none');
        gprMdl_dof5o = fitrgp(x5,y5,'Fitmethod','exact','PredictMethod','exact','Sigma',0.06,'KernelParameters',[0.15;0.35],'Basisfunction','none');
        gprMdl_dof6o = fitrgp(x6,y6,'Fitmethod','exact','PredictMethod','exact','Sigma',0.06,'KernelParameters',[0.15;0.35],'Basisfunction','none');
end

% Hyperparameter matrix on the Gaussian process regression
hp = [gprMdl_dof1o.Sigma,gprMdl_dof1o.KernelInformation.KernelParameters'; ...
      gprMdl_dof2o.Sigma,gprMdl_dof2o.KernelInformation.KernelParameters'; ...
      gprMdl_dof3o.Sigma,gprMdl_dof3o.KernelInformation.KernelParameters'; ...
      gprMdl_dof4o.Sigma,gprMdl_dof4o.KernelInformation.KernelParameters'; ...
      gprMdl_dof5o.Sigma,gprMdl_dof5o.KernelInformation.KernelParameters'; ...
      gprMdl_dof6o.Sigma,gprMdl_dof6o.KernelInformation.KernelParameters'];
  
% Make hyper parameters to build GPR model. Baumeister chose these based on tests.
hp_new = [hp(:,1)*(0.1/mean(hp(:,1))),hp(:,2)*(0.7/mean(hp(:,2))),hp(:,3)*(0.9/mean(hp(:,3)))];

% GPR with adjusted kernel parameters. 
gprMdl_dof1 = fitrgp(x1,y1,'Fitmethod','none','Sigma',0.1,'KernelParameters',[hp_new(1,2);hp_new(1,3)],'Basisfunction','none');
gprMdl_dof2 = fitrgp(x2,y2,'Fitmethod','none','Sigma',0.1,'KernelParameters',[hp_new(2,2);hp_new(2,3)],'Basisfunction','none');
gprMdl_dof3 = fitrgp(x3,y3,'Fitmethod','none','Sigma',0.1,'KernelParameters',[hp_new(3,2);hp_new(3,3)],'Basisfunction','none');
gprMdl_dof4 = fitrgp(x4,y4,'Fitmethod','none','Sigma',0.1,'KernelParameters',[hp_new(4,2);hp_new(4,3)],'Basisfunction','none');
gprMdl_dof5 = fitrgp(x5,y5,'Fitmethod','none','Sigma',0.1,'KernelParameters',[hp_new(4,2);hp_new(4,3)],'Basisfunction','none');
gprMdl_dof6 = fitrgp(x6,y6,'Fitmethod','none','Sigma',0.1,'KernelParameters',[hp_new(4,2);hp_new(4,3)],'Basisfunction','none');

%% LR
% LR - basic model (intersection incl)
LRmdl_1 = fitlm(x1,y1);
LRmdl_2 = fitlm(x2,y2);
LRmdl_3 = fitlm(x3,y3);
LRmdl_4 = fitlm(x4,y4);
LRmdl_5 = fitlm(x5,y5);
LRmdl_6 = fitlm(x6,y6);

%% Test predictions on training data
% GPR
% the training data
[y1_testGPR,p1,o1] = predict(gprMdl_dof1,RMS);
[y2_testGPR,p2,o2] = predict(gprMdl_dof2,RMS);
[y3_testGPR,p3,o3] = predict(gprMdl_dof3,RMS);
[y4_testGPR,p4,o4] = predict(gprMdl_dof4,RMS);
[y5_testGPR,p5,o5] = predict(gprMdl_dof5,RMS);
[y6_testGPR,p6,o6] = predict(gprMdl_dof6,RMS);

%% Plots
% GPR
subplot = @(m,n,p) subtightplot(m, n, p, [0.06 0.03], [0.08 0.03], [0.08 0.03]);
LR_GPR=figure('DefaultAxesPosition', [0.1, 0.1, 0.8, 0.8], 'units', 'normalized','position', [0 0 0.39774756441 1],'Color','w');
ax1 = subplot(4,1,1);
rectangle('Position',[i(3) -0.24 i(2)-1 1.25],'FaceColor',[.95,.95,.95],'Linestyle','none');
rectangle('Position',[0 -0.25 i(2)-1 1.25],'EdgeColor',[.4,.4,.4],'Linewidth',1.2); hold on;
plot(1:length(RMS),y1_testGPR,'k');
title('Flexion')
hold on;
ax2 = subplot(4,1,2);
rectangle('Position',[i(3) -0.23 i(2)-1 1.25],'FaceColor',[.95,.95,.95],'Linestyle','none');
rectangle('Position',[i(2) -0.25 i(2)-1 1.25],'EdgeColor',[.4,.4,.4],'Linewidth',1.2); hold on;
plot(1:length(RMS),y2_testGPR,'k');
title('Extension')
ylabel('normalized [-]')
ax3 = subplot(4,1,3);
rectangle('Position',[i(3) -0.23 i(2)-1 1.25],'FaceColor',[.95,.95,.95],'Linestyle','none');
rectangle('Position',[i(4) -0.25 i(2)-1 1.25],'EdgeColor',[.4,.4,.4],'Linewidth',1.2); hold on;
plot(1:length(RMS),y3_testGPR,'k');
title('Abduction')
ax4 = subplot(4,1,4);
rectangle('Position',[i(3) -0.23 i(2)-1 1.25],'FaceColor',[.95,.95,.95],'Linestyle','none');
rectangle('Position',[i(5) -0.25 i(2)-1 1.25],'EdgeColor',[.4,.4,.4],'Linewidth',1.2); hold on;
plot(1:length(RMS),y4_testGPR,'k');
title('Adduction')
xlabel('samples [#]')
[~, hobj, ~, ~] = legend('LR','GPR','Location','Best');
set(hobj,'linewidth',1.5);
gca_handles = [ax1,ax2,ax3,ax4];
set(gca_handles,'fontsize',10,'YLim',[-0.25 1],'XLim',[0 length(y1_testGPR)+1])

% P-values (std) of GPR
subplot = @(m,n,p) subtightplot(m, n, p, [0.06 0.03], [0.08 0.05], [0.08 0.03]);
GPR_p=figure('DefaultAxesPosition', [0.1, 0.1, 0.8, 0.8], 'units', 'normalized','position', [0 0 0.39774756441 1],'Color','w');
ax1=subplot(4,1,1);
rectangle('Position',[i(3) 0.008 i(2)-1 0.39],'FaceColor',[.95,.95,.95],'Linestyle','none');
rectangle('Position',[0 0 i(2)-1 0.1],'EdgeColor',[.4,.4,.4],'Linewidth',1.2); hold on;
plot(1:length(RMS),p1);
title('Flexion')
hold on;
ax2=subplot(4,1,2);
rectangle('Position',[i(3) 0.008 i(2)-1 0.39],'FaceColor',[.95,.95,.95],'Linestyle','none');
rectangle('Position',[i(2) 0 i(2)-1 0.1],'EdgeColor',[.4,.4,.4],'Linewidth',1.2); hold on;
plot(1:length(RMS),p2);
title('Extension')
ax3=subplot(4,1,3);
rectangle('Position',[i(3) 0.008 i(2)-1 0.39],'FaceColor',[.95,.95,.95],'Linestyle','none');
rectangle('Position',[i(4) 0 i(2)-1 0.1],'EdgeColor',[.4,.4,.4],'Linewidth',1.2); hold on;
plot(1:length(RMS),p3);
title('Abduction')
ylabel('uncertainty [ratio]')
ax4=subplot(4,1,4);
rectangle('Position',[i(3) 0.008 i(2)-1 0.39],'FaceColor',[.95,.95,.95],'Linestyle','none');
rectangle('Position',[i(5) 0 i(2)-3 0.1],'EdgeColor',[.4,.4,.4],'Linewidth',1.2); hold on;
plot(1:length(RMS),p4);
title('Adduction')
xlabel('samples [#]')
[~, hobj, ~, ~] = legend('p_{GPR}','Location','Best');
set(hobj,'linewidth',1.5);
gca_handles = [ax1,ax2,ax3,ax4];
set(gca_handles,'fontsize',10,'YLim',[0 0.1],'XLim',[0 length(y1_testGPR)+1])

%% Save Figures
