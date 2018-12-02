clc; clear; close all;
%% ST7 - Gaussian Process Regression (GPR)
%This script uses data collected using the MATLAB-app
%'CollectTrainingData', to build a GPR-model.

%% Loading and structuring data
% Prompt the user to choose the folder containing the data of the
% respective test subject.
fileFolder = dir(uigetdir);
fileFolder = fileFolder.folder;
addpath(fileFolder);

% Get a list of all .mat-files in the folder.
listOfFiles = dir(fullfile(fileFolder, '*.mat'));
% Load the files
EMG= []; MVC=[]; RMS= []; GeneratedProfile = []; movement=[];

%% Choose number of degrees of freedom
        
        for k = 1 : length(listOfFiles)
            
            load(listOfFiles(k).name, 'emg');
            data{1,k} = emg.values;
            data{2,k} = emg.MVC;
            EMG = [EMG; data{1,k}];
            MVC = [MVC; data{2,k}];
            
            RMS_temp = [];
            for i=1 : 8
                RMS_temp(:,i) = transpose(my_rms(data{1,k}(:,i),40,20,0));
            end
            
            RMS = [RMS; RMS_temp];
            
            GeneratedProfile = [GeneratedProfile my_rms(mean(data{1,k},2), 40, 20, 0)];
            
            %To be used for indexing - the lengths of the RMS from each movement:
            f(k) = length(RMS_temp);
            
        end
        clear emg RMS_temp;
        
        GeneratedProfile= transpose(GeneratedProfile);
        
        %Assigning a number for each movement to be used in indexing. Feasible
        %because the files are always loaded alphabetically.
        movement(1:sum(f(1:3))) =2;                   %Close = 1
        movement(sum(f(1:3))+1: sum(f(1:6)))=3;       %Extension= 2
        movement(sum(f(1:6))+1:sum(f(1:9)))=5;        %Flexion = 3
        movement(sum(f(1:9))+1:sum(f(1:12)))=6;       %Open = 4
        movement(sum(f(1:12))+1:sum(f(1:15)))=7;      %Radial deviation = 5
        
        %The index vector
        [dofs, i, ~] = unique(movement);
        MVC = unique(MVC, 'stable');
        
        fig_rawData= figure('DefaultAxesPosition', [0.1, 0.1, 0.8, 0.8], 'units', 'normalized','position', [0 0 0.8 1],'Color','w');
        subplot(3,1,1)
        plot(1:length(EMG),EMG)
        subplot(3,1,2)
        plot(1:length(RMS),RMS)
        subplot(3,1,3)
        plot(1:length(GeneratedProfile), GeneratedProfile)
        export_fig(fig_rawData,'-pdf','filename','plot_rawData');
        %% Arrangement in predictors and target
        
        %Predictor - RMS of each EMG-channel of the myo armband for 25%, 50%, 75% of MVC.
        %Vi tr?kker 1 fra fordi i(2) er f?rste sample i n?ste bev?gelse.
        x_ext = RMS(i(1):i(2)-1,:);
        x_flex = RMS(i(2):i(3)-1,:);
        x_rd = RMS(i(3):i(4)-1,:);
        x_rest = RMS(i(4):i(5)-1,:);
        x_ud = RMS(i(5):end,:);
        %x_ud = RMS(i(7):i(8)-1,:);
        
        % Setting rest to zero
        GeneratedProfile(i(4):i(5)-1,:)= 0;
        %GeneratedProfile(i(8):end,:)= 0;
        
        %Target values - The generated profile.
        y_ext = GeneratedProfile(i(1):i(2)-1,:);
        y_flex = GeneratedProfile(i(2):i(3)-1,:);
        y_rd = GeneratedProfile(i(3):i(4)-1,:);
        y_rest = GeneratedProfile(i(4):i(5)-1,:);
        y_ud = GeneratedProfile(i(5):end,:);
        %y_ud = GeneratedProfile(i(7):i(8)-1,:);
        
 %Opsplitning til brug i træning 'andre nul'
                y_reg = [y_ext,zeros(length(y_ext),4);...
                    zeros(length(y_flex),1), y_flex, zeros(length(y_flex),3);
                    zeros(length(y_rd),2),y_rd, zeros(length(y_rd),2);
                    zeros(length(y_rest),3),y_rest, zeros(length(y_rest),1);
                    zeros(length(y_ud),4),y_ud];
                
                y1_oz = y_reg(:,1); %ext oz=others zero
                y2_oz = y_reg(:,2); %flex
                y3_oz = y_reg(:,3); %rd
                y4_oz = y_reg(:,4); %rest
                y5_oz= y_reg(:,5); %ud
                x1_oz = RMS; x2_oz = RMS; x3_oz = RMS; x4_oz = RMS; x5_oz =RMS;
                
 %Opsplitning til brug i træning 'individuelt'.
                y1_indi = [y_ext;y_rest];    x1_indi = [x_ext;x_rest];
                y2_indi = [y_flex;y_rest];   x2_indi = [x_flex;x_rest];
                y3_indi = [y_rd;y_rest];     x3_indi = [x_rd;x_rest];
                y5_indi = [y_ud;y_rest];     x5_indi = [x_ud;x_rest];
     
 %% Gaussian Process Regression
        
%GPR trænet med andre til nul (oz) til brug i GPR u. confidence
gprMdl_dof2_uc = fitrgp(x1_oz,y1_oz,'Fitmethod','sr','PredictMethod','sr','Sigma',0.06,'KernelParameters',[0.15;0.35],'Basisfunction','none');
gprMdl_dof3_uc = fitrgp(x2_oz,y2_oz,'Fitmethod','sr','PredictMethod','sr','Sigma',0.06,'KernelParameters',[0.15;0.35],'Basisfunction','none');
gprMdl_dof5_uc = fitrgp(x3_oz,y3_oz,'Fitmethod','sr','PredictMethod','sr','Sigma',0.06,'KernelParameters',[0.15;0.35],'Basisfunction','none');
gprMdl_dof6_uc = fitrgp(x5_oz,y5_oz,'Fitmethod','sr','PredictMethod','sr','Sigma',0.06,'KernelParameters',[0.15;0.35],'Basisfunction','none');
     
hp_oz = [gprMdl_dof2_uc.Sigma,gprMdl_dof2_uc.KernelInformation.KernelParameters'; ...
      gprMdl_dof3_uc.Sigma,gprMdl_dof3_uc.KernelInformation.KernelParameters'; ...
      gprMdl_dof5_uc.Sigma,gprMdl_dof5_uc.KernelInformation.KernelParameters'; ...
      gprMdl_dof6_uc.Sigma,gprMdl_dof6_uc.KernelInformation.KernelParameters'];      
        
% Make hyper parameters to build GPR model. Chosen based on tests.
hp_new_oz = [hp_oz(:,1)*(0.1/mean(hp_oz(:,1))),hp_oz(:,2)*(0.7/mean(hp_oz(:,2))),hp_oz(:,3)*(0.7/mean(hp_oz(:,3)))];

% GPR with adjusted kernel parameters. This one is used for the experiment
gprMdl_dof2_uc = fitrgp(x1_oz,y1_oz,'Fitmethod','none','Sigma',0.01,'KernelParameters',[hp_new_oz(1,2);hp_new_oz(1,3)],'Basisfunction','none');
gprMdl_dof3_uc = fitrgp(x2_oz,y2_oz,'Fitmethod','none','Sigma',0.01,'KernelParameters',[hp_new_oz(2,2);hp_new_oz(2,3)],'Basisfunction','none');
gprMdl_dof5_uc = fitrgp(x3_oz,y3_oz,'Fitmethod','none','Sigma',0.01,'KernelParameters',[hp_new_oz(3,2);hp_new_oz(3,3)],'Basisfunction','none');
gprMdl_dof6_uc = fitrgp(x5_oz,y5_oz,'Fitmethod','none','Sigma',0.01,'KernelParameters',[hp_new_oz(4,2);hp_new_oz(4,3)],'Basisfunction','none');


%GPR trænet individuelt (indi) til brug i GPR m. confidence
gprMdl_dof2_mc = fitrgp(x1_indi,y1_indi,'Fitmethod','sr','PredictMethod','sr','Sigma',0.06,'KernelParameters',[0.15;0.35],'Basisfunction','none');
gprMdl_dof3_mc = fitrgp(x2_indi,y2_indi,'Fitmethod','sr','PredictMethod','sr','Sigma',0.06,'KernelParameters',[0.15;0.35],'Basisfunction','none');
gprMdl_dof5_mc = fitrgp(x3_indi,y3_indi,'Fitmethod','sr','PredictMethod','sr','Sigma',0.06,'KernelParameters',[0.15;0.35],'Basisfunction','none');
gprMdl_dof6_mc = fitrgp(x5_indi,y5_indi,'Fitmethod','sr','PredictMethod','sr','Sigma',0.06,'KernelParameters',[0.15;0.35],'Basisfunction','none');
     
hp_indi = [gprMdl_dof2_mc.Sigma,gprMdl_dof2_mc.KernelInformation.KernelParameters'; ...
      gprMdl_dof3_mc.Sigma,gprMdl_dof3_mc.KernelInformation.KernelParameters'; ...
      gprMdl_dof5_mc.Sigma,gprMdl_dof5_mc.KernelInformation.KernelParameters'; ...
      gprMdl_dof6_mc.Sigma,gprMdl_dof6_mc.KernelInformation.KernelParameters'];      
        
% Make hyper parameters to build GPR model. Chosen based on tests.
hp_new_indi = [hp_indi(:,1)*(0.1/mean(hp_indi(:,1))),hp_indi(:,2)*(0.7/mean(hp_indi(:,2))),hp_indi(:,3)*(0.7/mean(hp_indi(:,3)))];

% GPR with adjusted kernel parameters. This one is used for the experiment
gprMdl_dof2_mc = fitrgp(x1_indi,y1_indi,'Fitmethod','none','Sigma',0.01,'KernelParameters',[hp_new_indi(1,2);hp_new_indi(1,3)],'Basisfunction','none');
gprMdl_dof3_mc = fitrgp(x2_indi,y2_indi,'Fitmethod','none','Sigma',0.01,'KernelParameters',[hp_new_indi(2,2);hp_new_indi(2,3)],'Basisfunction','none');
gprMdl_dof5_mc = fitrgp(x3_indi,y3_indi,'Fitmethod','none','Sigma',0.01,'KernelParameters',[hp_new_indi(3,2);hp_new_indi(3,3)],'Basisfunction','none');
gprMdl_dof6_mc = fitrgp(x5_indi,y5_indi,'Fitmethod','none','Sigma',0.01,'KernelParameters',[hp_new_indi(4,2);hp_new_indi(4,3)],'Basisfunction','none');


%% LR
LRmdl_2 = fitlm(x1_oz,y1_oz);
LRmdl_3 = fitlm(x2_oz,y2_oz);
LRmdl_5 = fitlm(x3_oz,y3_oz);
LRmdl_6 = fitlm(x5_oz,y5_oz);
   