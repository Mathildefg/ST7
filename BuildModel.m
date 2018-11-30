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
train_dof = input('#DoF? (2 or 3):');
train_type = input('What training? Set others zero (1) or Individual (2) Opposite zero (3): ');
switch train_dof
    %% Dof 2
    case 2
        
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
        %because the files are always loaded alphabetically9
        movement(1:sum(f(1:3))) =1;                   %Close = 1
        movement(sum(f(1:3))+1: sum(f(1:6)))=2;       %Extension= 2
        movement(sum(f(1:6))+1:sum(f(1:9)))=3;        %Flexion = 3
        movement(sum(f(1:9))+1:sum(f(1:12)))=4;       %Open = 4
        movement(sum(f(1:12))+1:sum(f(1:15)))=5;      %Radial deviation = 5
        movement(sum(f(1:15))+1:sum(f(1:18)))= 6;     %Rest = 6
        movement(sum(f(1:18))+1:sum(f(1:21)))=7;      %Ulnar deviation = 7
        
        %The index vector
        [dofs, i, ~] = unique(movement);
        MVC = unique(MVC, 'stable');
        
        figure(1)
        subplot(3,1,1)
        plot(1:length(EMG),EMG)
        subplot(3,1,2)
        plot(1:length(RMS),RMS)
        subplot(3,1,3)
        plot(1:length(GeneratedProfile), GeneratedProfile)
        %% Arrangement in predictors and target
        
        %Predictor - RMS of each EMG-channel of the myo armband for 25%, 50%, 75% of MVC.
        x_cls = RMS(i(1):i(2)-1,:);  %Vi tr?kker 1 fra fordi i(2) er f?rste sample i n?ste bev?gelse.
        x_ext = RMS(i(2):i(3)-1,:);
        x_flex = RMS(i(3):i(4)-1,:);
        x_opn = RMS(i(4):i(5)-1,:);
        x_rd = RMS(i(5):i(6)-1,:);
        x_rest = RMS(i(6):i(7)-1,:);
        x_ud = RMS(i(7):end,:);
        %x_ud = RMS(i(7):i(8)-1,:);
        
        % Setting rest to zero
        GeneratedProfile(i(6):i(7)-1,:)= 0;
        %GeneratedProfile(i(8):end,:)= 0;
        
        %Target values - The generated profile.
        y_cls = GeneratedProfile(i(1):i(2)-1,:);
        y_ext = GeneratedProfile(i(2):i(3)-1,:);
        y_flex = GeneratedProfile(i(3):i(4)-1,:);
        y_opn = GeneratedProfile(i(4):i(5)-1,:);
        y_rd = GeneratedProfile(i(5):i(6)-1,:);
        y_rest = GeneratedProfile(i(6):i(7)-1,:);
        y_ud = GeneratedProfile(i(7):end,:);
        %y_ud = GeneratedProfile(i(7):i(8)-1,:);
        switch train_type
            case 1
                y_reg = [y_ext,zeros(length(y_ext),4);...
                    zeros(length(y_flex),1), y_flex, zeros(length(y_flex),3);
                    zeros(length(y_rd),2),y_rd, zeros(length(y_rd),2);
                    zeros(length(y_rest),3),y_rest, zeros(length(y_rest),1);
                    zeros(length(y_ud),4),y_ud];
                
                y1 = y_reg(:,1); %ext
                y2 = y_reg(:,2); %flex
                y3 = y_reg(:,3); %rd
                y4 = y_reg(:,4); %rest
                y5= y_reg(:,5); %ud
                x1 = RMS2d; x2 = RMS2d; x3 = RMS2d; x4 = RMS2d; x5=RMS2d;
            case 2
                y1 = [y_ext;y_rest];    x1 = [x_ext;x_rest];
                y2 = [y_flex;y_rest];   x2 = [x_flex;x_rest];
                y3 = [y_rd;y_rest];     x3 = [x_rd;x_rest];
                y5 = [y_ud;y_rest];     x5 = [x_ud;x_rest];
            case 3
                y1= [y_ext;zeros(length(y_flex),1);zeros(length(y_rest),1)];
                y2= [y_flex;zeros(length(y_ext),1);zeros(length(y_rest),1)];
                y3= [y_rd;zeros(length(y_ud),1);zeros(length(y_rest),1)];
                y5= [y_ud;zeros(length(y_rd),1);zeros(length(y_rest),1)];
                
                x1 = [x_ext;x_flex;x_rest];
                x2 = [x_flex;x_ext;x_rest];
                x3 = [x_rd;x_ud;x_rest];
                x5 = [x_ud;x_rd;x_rest];
        end
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
                gprMdl_dof2 = fitrgp(x1,y1,'FitMethod','exact','PredictMethod','exact','Basisfunction','none','OptimizeHyperparameters',{'Sigma'},'HyperparameterOptimizationOptions', struct('MaxObjectiveEvaluations',5));
                gprMdl_dof3 = fitrgp(x2,y2,'FitMethod','exact','PredictMethod','exact','Basisfunction','none','OptimizeHyperparameters',{'Sigma'},'HyperparameterOptimizationOptions', struct('MaxObjectiveEvaluations',5));
                gprMdl_dof5 = fitrgp(x3,y3,'FitMethod','exact','PredictMethod','exact','Basisfunction','none','OptimizeHyperparameters',{'Sigma'},'HyperparameterOptimizationOptions', struct('MaxObjectiveEvaluations',5));
                gprMdl_dof6 = fitrgp(x5,y5,'FitMethod','exact','PredictMethod','exact','Basisfunction','none','OptimizeHyperparameters',{'Sigma'},'HyperparameterOptimizationOptions', struct('MaxObjectiveEvaluations',5));
                
            case 'N'        % if non-optimization (already known hyperparameters) method is chosen
                % Gaussian Regression - Already known hyperparameters
                gprMdl_dof2 = fitrgp(x1,y1,'Fitmethod','sr','PredictMethod','sr','Sigma',0.06,'KernelParameters',[0.15;0.35],'Basisfunction','none');
                gprMdl_dof3 = fitrgp(x2,y2,'Fitmethod','sr','PredictMethod','sr','Sigma',0.06,'KernelParameters',[0.15;0.35],'Basisfunction','none');
                gprMdl_dof5 = fitrgp(x3,y3,'Fitmethod','sr','PredictMethod','sr','Sigma',0.06,'KernelParameters',[0.15;0.35],'Basisfunction','none');
                gprMdl_dof6 = fitrgp(x5,y5,'Fitmethod','sr','PredictMethod','sr','Sigma',0.06,'KernelParameters',[0.15;0.35],'Basisfunction','none');
        end
        
        
        
        %% LR
        % LR - basic model (intersection incl)
        LRmdl_2 = fitlm(x1,y1);
        LRmdl_3 = fitlm(x2,y2);
        LRmdl_5 = fitlm(x3,y3);
        LRmdl_6 = fitlm(x5,y5);
        
        %% Test predictions on training data
        % GPR
        % the training data
        [y2_testGPR,p2,o2] = predict(gprMdl_dof2,RMS2d);
        [y3_testGPR,p3,o3] = predict(gprMdl_dof3,RMS2d);
        [y5_testGPR,p5,o5] = predict(gprMdl_dof5,RMS2d);
        [y6_testGPR,p6,o6] = predict(gprMdl_dof6,RMS2d);
        
        % LR
        y2_LR = predict(LRmdl_2,RMS2d);
        y3_LR = predict(LRmdl_3,RMS2d);
        y5_LR = predict(LRmdl_5,RMS2d);
        y6_LR = predict(LRmdl_6,RMS2d);
        
        %% Plots
        % GPR
        subplot = @(m,n,p) subtightplot(m, n, p, [0.06 0.03], [0.08 0.03], [0.08 0.03]);
        LR_GPR=figure('DefaultAxesPosition', [0.1, 0.1, 0.8, 0.8], 'units', 'normalized','position', [0 0 0.39774756441 1],'Color','w');
        ax1 = subplot(6,1,1);
        rectangle('Position',[i(4) -0.25 i(2)-1 1.25],'FaceColor',[.95,.95,.95],'Linestyle','none');
        rectangle('Position',[0 -0.25 i(2)-1 1.25],'EdgeColor',[.4,.4,.4],'Linewidth',1.2); hold on;
        plot(1:length(RMS2d),y2_testGPR,'k');
        title('Extension')
        hold on;
        plot(1:length(RMS2d),y2_LR,'b');
        plot(1:length(RMS2d),o2,'--k');
        ax2 = subplot(6,1,2);
        rectangle('Position',[i(4) -0.25 i(2)-1 1.25],'FaceColor',[.95,.95,.95],'Linestyle','none');
        rectangle('Position',[i(2) -0.25 i(2)-1 1.25],'EdgeColor',[.4,.4,.4],'Linewidth',1.2); hold on;
        plot(1:length(RMS2d),y3_testGPR,'k');
        plot(1:length(RMS2d),y3_LR,'b');
        plot(1:length(RMS2d),o3,'--k');
        title('Flexion')
        ax3 = subplot(6,1,3);
        rectangle('Position',[i(4) -0.25 i(2)-1 1.25],'FaceColor',[.95,.95,.95],'Linestyle','none');
        rectangle('Position',[i(3) -0.25 i(2)-1 1.25],'EdgeColor',[.4,.4,.4],'Linewidth',1.2); hold on;
        plot(1:length(RMS2d),y5_testGPR,'k');
        plot(1:length(RMS2d),y5_LR,'b');
        plot(1:length(RMS2d),o5,'--k');
        title('Radial deviation')
        ylabel('normalized [-]')
        ax4 = subplot(6,1,4);
        rectangle('Position',[i(4) -0.25 i(2)-1 1.25],'FaceColor',[.95,.95,.95],'Linestyle','none');
        rectangle('Position',[i(5) -0.23 i(2)-1 1.25],'EdgeColor',[.4,.4,.4],'Linewidth',1.2); hold on;
        plot(1:length(RMS2d),y6_testGPR,'k');
        plot(1:length(RMS2d),y6_LR,'b');
        plot(1:length(RMS2d),o6,'--k');
        title('Open')
        title('Ulna deviation')
        xlabel('samples [#]')
        [~, hobj, ~, ~] = legend('GPR','LR','CI(GPR)','Location','Best');
        set(hobj,'linewidth',1.5);
        gca_handles = [ax1,ax2,ax3,ax4];
        set(gca_handles,'fontsize',10,'YLim',[-0.25 0.5],'XLim',[0 length(y2_testGPR)+1])
    case 3
        %% Dof 3
        
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
        %because the files are always loaded alphabetically9
        movement(1:sum(f(1:3))) =1;                   %Close = 1
        movement(sum(f(1:3))+1: sum(f(1:6)))=2;       %Extension= 2
        movement(sum(f(1:6))+1:sum(f(1:9)))=3;        %Flexion = 3
        movement(sum(f(1:9))+1:sum(f(1:12)))=4;       %Open = 4
        movement(sum(f(1:12))+1:sum(f(1:15)))=5;      %Radial deviation = 5
        movement(sum(f(1:15))+1:sum(f(1:18)))= 6;     %Rest = 6
        movement(sum(f(1:18))+1:sum(f(1:21)))=7;      %Ulnar deviation = 7
        
        %The index vector
        [dofs, i, ~] = unique(movement);
        MVC = unique(MVC, 'stable');
        
        figure(1)
        subplot(3,1,1)
        plot(1:length(EMG),EMG)
        subplot(3,1,2)
        plot(1:length(RMS),RMS)
        subplot(3,1,3)
        plot(1:length(GeneratedProfile), GeneratedProfile)
        %% Arrangement in predictors and target
        
        %Predictor - RMS of each EMG-channel of the myo armband for 25%, 50%, 75% of MVC.
        x_cls = RMS(i(1):i(2)-1,:);  %Vi tr?kker 1 fra fordi i(2) er f?rste sample i n?ste bev?gelse.
        x_ext = RMS(i(2):i(3)-1,:);
        x_flex = RMS(i(3):i(4)-1,:);
        x_opn = RMS(i(4):i(5)-1,:);
        x_rd = RMS(i(5):i(6)-1,:);
        x_rest = RMS(i(6):i(7)-1,:);
        x_ud = RMS(i(7):end,:);
        %x_ud = RMS(i(7):i(8)-1,:);
        
        % Setting rest to zero
        GeneratedProfile(i(6):i(7)-1,:)= 0;
        %GeneratedProfile(i(8):end,:)= 0;
        
        %Target values - The generated profile.
        y_cls = GeneratedProfile(i(1):i(2)-1,:);
        y_ext = GeneratedProfile(i(2):i(3)-1,:);
        y_flex = GeneratedProfile(i(3):i(4)-1,:);
        y_opn = GeneratedProfile(i(4):i(5)-1,:);
        y_rd = GeneratedProfile(i(5):i(6)-1,:);
        y_rest = GeneratedProfile(i(6):i(7)-1,:);
        y_ud = GeneratedProfile(i(7):end,:);
        %y_ud = GeneratedProfile(i(7):i(8)-1,:);
        %% Choose between training individual or with other setting to zero
        % Setting others to 0 is performed during the experiment and has shown to
        % derive better performance by tests.
        % if train_type ~= 2
        %     train_type = 1;
        % end
        
        switch train_type
            case 1
                y_reg = [y_cls,zeros(length(y_cls),6);
                    zeros(length(y_ext),1),y_ext,zeros(length(y_ext),5);...
                    zeros(length(y_flex),2), y_flex, zeros(length(y_flex),4);
                    zeros(length(y_opn),3),y_opn, zeros(length(y_opn),3);...
                    zeros(length(y_rd),4),y_rd, zeros(length(y_rd),2);
                    zeros(length(y_rest),5),y_rest, zeros(length(y_rest),1);
                    zeros(length(y_ud),6),y_ud];
                
                y1 = y_reg(:,1); %cls
                y2 = y_reg(:,2); %ext
                y3 = y_reg(:,3); %flex
                y4 = y_reg(:,4); %opn
                y5= y_reg(:,5); %rd
                y6= y_reg(:,7); %ud
                x1 = RMS; x2 = RMS; x3 = RMS; x4 = RMS; x5=RMS; x6=RMS;
                
            case 2
                y1 = [y_cls;y_rest];   x1 = [x_cls;x_rest];
                y2 = [y_ext;y_rest];   x2 = [x_ext;x_rest];
                y3 = [y_flex;y_rest];  x3 = [x_flex;x_rest];
                y4 = [y_opn;y_rest];   x4 = [x_opn;x_rest];
                y5 = [y_rd;y_rest];    x5 = [x_rd;x_rest];
                y6 = [y_ud;y_rest];    x6 = [x_ud;x_rest];
                
                
            case 3 %Opposite zero
                
                y1= [y_cls;zeros(length(y_opn),1);zeros(length(y_rest),1)];
                y2= [y_ext;zeros(length(y_flex),1);zeros(length(y_rest),1)];
                y3= [y_flex;zeros(length(y_ext),1);zeros(length(y_rest),1)];
                y4= [y_opn;zeros(length(y_cls),1);zeros(length(y_rest),1)];
                y5= [y_rd;zeros(length(y_ud),1);zeros(length(y_rest),1)];
                y6= [y_ud;zeros(length(y_rd),1);zeros(length(y_rest),1)];
                
                x1 = [x_cls;x_opn;x_rest];
                x2 = [x_ext;x_flex;x_rest];
                x3 = [x_flex;x_ext;x_rest];
                x4 = [x_opn;x_cls;x_rest];
                x5 = [x_rd;x_ud;x_rest];
                x6 = [x_ud;x_rd;x_rest];
                
        end
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
                gprMdl_dof1o = fitrgp(x1,y1,'FitMethod','exact','PredictMethod','exact','Basisfunction','none','OptimizeHyperparameters',{'Sigma'},'HyperparameterOptimizationOptions', struct('MaxObjectiveEvaluations',5));
                gprMdl_dof2o = fitrgp(x2,y2,'FitMethod','exact','PredictMethod','exact','Basisfunction','none','OptimizeHyperparameters',{'Sigma'},'HyperparameterOptimizationOptions', struct('MaxObjectiveEvaluations',5));
                gprMdl_dof3o = fitrgp(x3,y3,'FitMethod','exact','PredictMethod','exact','Basisfunction','none','OptimizeHyperparameters',{'Sigma'},'HyperparameterOptimizationOptions', struct('MaxObjectiveEvaluations',5));
                gprMdl_dof4o = fitrgp(x4,y4,'FitMethod','exact','PredictMethod','exact','Basisfunction','none','OptimizeHyperparameters',{'Sigma'},'HyperparameterOptimizationOptions', struct('MaxObjectiveEvaluations',5));
                gprMdl_dof5o = fitrgp(x5,y5,'FitMethod','exact','PredictMethod','exact','Basisfunction','none','OptimizeHyperparameters',{'Sigma'},'HyperparameterOptimizationOptions', struct('MaxObjectiveEvaluations',5));
                gprMdl_dof6o = fitrgp(x6,y6,'FitMethod','exact','PredictMethod','exact','Basisfunction','none','OptimizeHyperparameters',{'Sigma'},'HyperparameterOptimizationOptions', struct('MaxObjectiveEvaluations',5));
                
            case 'N'        % if non-optimization (already known hyperparameters) method is chosen
                % Gaussian Regression - Already known hyperparameters
                gprMdl_dof1o = fitrgp(x1,y1,'Fitmethod','sr','PredictMethod','sr','Sigma',0.06,'KernelParameters',[0.15;0.35],'Basisfunction','none');
                gprMdl_dof2o = fitrgp(x2,y2,'Fitmethod','sr','PredictMethod','sr','Sigma',0.06,'KernelParameters',[0.15;0.35],'Basisfunction','none');
                gprMdl_dof3o = fitrgp(x3,y3,'Fitmethod','sr','PredictMethod','sr','Sigma',0.06,'KernelParameters',[0.15;0.35],'Basisfunction','none');
                gprMdl_dof4o = fitrgp(x4,y4,'Fitmethod','sr','PredictMethod','sr','Sigma',0.06,'KernelParameters',[0.15;0.35],'Basisfunction','none');
                gprMdl_dof5o = fitrgp(x5,y5,'Fitmethod','sr','PredictMethod','sr','Sigma',0.06,'KernelParameters',[0.15;0.35],'Basisfunction','none');
                gprMdl_dof6o = fitrgp(x6,y6,'Fitmethod','sr','PredictMethod','sr','Sigma',0.06,'KernelParameters',[0.15;0.35],'Basisfunction','none');
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
        
        %GPR with adjusted kernel parameters.
        gprMdl_dof1 = fitrgp(x1,y1,'Fitmethod','none','Sigma',0.1,'KernelParameters',[hp_new(1,2);hp_new(1,3)],'Basisfunction','none');
        gprMdl_dof2 = fitrgp(x2,y2,'Fitmethod','none','Sigma',0.1,'KernelParameters',[hp_new(2,2);hp_new(2,3)],'Basisfunction','none');
        gprMdl_dof3 = fitrgp(x3,y3,'Fitmethod','none','Sigma',0.1,'KernelParameters',[hp_new(3,2);hp_new(3,3)],'Basisfunction','none');
        gprMdl_dof4 = fitrgp(x4,y4,'Fitmethod','none','Sigma',0.1,'KernelParameters',[hp_new(4,2);hp_new(4,3)],'Basisfunction','none');
        gprMdl_dof5 = fitrgp(x5,y5,'Fitmethod','none','Sigma',0.1,'KernelParameters',[hp_new(5,2);hp_new(5,3)],'Basisfunction','none');
        gprMdl_dof6 = fitrgp(x6,y6,'Fitmethod','none','Sigma',0.1,'KernelParameters',[hp_new(6,2);hp_new(6,3)],'Basisfunction','none');
        
        % GPR model with fixed hyper parameters.
        %  gprMdl_dof1 = fitrgp(x1,y1,'Fitmethod','none','Sigma',0.1,'KernelParameters',[0.7;0.9],'Basisfunction','none');
        %  gprMdl_dof2 = fitrgp(x2,y2,'Fitmethod','none','Sigma',0.1,'KernelParameters',[0.7;0.9],'Basisfunction','none');
        %  gprMdl_dof3 = fitrgp(x3,y3,'Fitmethod','none','Sigma',0.1,'KernelParameters',[0.7;0.9],'Basisfunction','none');
        %  gprMdl_dof4 = fitrgp(x4,y4,'Fitmethod','none','Sigma',0.1,'KernelParameters',[0.7;0.9],'Basisfunction','none');
        %  gprMdl_dof5 = fitrgp(x5,y5,'Fitmethod','none','Sigma',0.1,'KernelParameters',[0.7;0.9],'Basisfunction','none');
        %  gprMdl_dof6 = fitrgp(x6,y6,'Fitmethod','none','Sigma',0.1,'KernelParameters',[0.7;0.9],'Basisfunction','none');
        
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
        [y1_testGPR,p1,o1] = predict(gprMdl_dof1,RMS);  %[predicted response values, estimated standard deviation, prediction intervals]
        [y2_testGPR,p2,o2] = predict(gprMdl_dof2,RMS);
        [y3_testGPR,p3,o3] = predict(gprMdl_dof3,RMS);
        [y4_testGPR,p4,o4] = predict(gprMdl_dof4,RMS);
        [y5_testGPR,p5,o5] = predict(gprMdl_dof5,RMS);
        [y6_testGPR,p6,o6] = predict(gprMdl_dof6,RMS);
        
        % LR
        y1_LR = predict(LRmdl_1,RMS);
        y2_LR = predict(LRmdl_2,RMS);
        y3_LR = predict(LRmdl_3,RMS);
        y4_LR = predict(LRmdl_4,RMS);
        y5_LR = predict(LRmdl_5,RMS);
        y6_LR = predict(LRmdl_6,RMS);
        
        %% Plots
        % GPR
        subplot = @(m,n,p) subtightplot(m, n, p, [0.06 0.03], [0.08 0.03], [0.08 0.03]);
        LR_GPR=figure('DefaultAxesPosition', [0.1, 0.1, 0.8, 0.8], 'units', 'normalized','position', [0 0 0.39774756441 1],'Color','w');
        ax1 = subplot(6,1,1);
        rectangle('Position',[i(6) -0.24 i(2)-1 1.25],'FaceColor',[.95,.95,.95],'Linestyle','none');
        rectangle('Position',[0 -0.25 i(2)-1 1.25],'EdgeColor',[.4,.4,.4],'Linewidth',1.2); hold on;
        plot(1:length(RMS),y1_testGPR,'k');
        title('Close')
        hold on;
        plot(1:length(RMS),y1_LR,'b');
        plot(1:length(RMS),o1,'--k');
        ax2 = subplot(6,1,2);
        rectangle('Position',[i(6) -0.23 i(2)-1 1.25],'FaceColor',[.95,.95,.95],'Linestyle','none');
        rectangle('Position',[i(2) -0.25 i(2)-1 1.25],'EdgeColor',[.4,.4,.4],'Linewidth',1.2); hold on;
        plot(1:length(RMS),y2_testGPR,'k');
        plot(1:length(RMS),y2_LR,'b');
        plot(1:length(RMS),o2,'--k');
        title('Extension')
        ax3 = subplot(6,1,3);
        rectangle('Position',[i(6) -0.23 i(2)-1 1.25],'FaceColor',[.95,.95,.95],'Linestyle','none');
        rectangle('Position',[i(3) -0.25 i(2)-1 1.25],'EdgeColor',[.4,.4,.4],'Linewidth',1.2); hold on;
        plot(1:length(RMS),y3_testGPR,'k');
        plot(1:length(RMS),y3_LR,'b');
        plot(1:length(RMS),o3,'--k');
        title('Flexion')
        ylabel('normalized [-]')
        ax4 = subplot(6,1,4);
        rectangle('Position',[i(6) -0.23 i(2)-1 1.25],'FaceColor',[.95,.95,.95],'Linestyle','none');
        rectangle('Position',[i(4) -0.25 i(2)-1 1.25],'EdgeColor',[.4,.4,.4],'Linewidth',1.2); hold on;
        plot(1:length(RMS),y4_testGPR,'k');
        plot(1:length(RMS),y4_LR,'b');
        plot(1:length(RMS),o4,'--k');
        title('Open')
        ax5 = subplot(6,1,5);
        rectangle('Position',[i(6) -0.23 i(2)-1 1.25],'FaceColor',[.95,.95,.95],'Linestyle','none');
        rectangle('Position',[i(5) -0.25 i(2)-1 1.25],'EdgeColor',[.4,.4,.4],'Linewidth',1.2); hold on;
        plot(1:length(RMS),y5_testGPR,'k');
        plot(1:length(RMS),y5_LR,'b');
        plot(1:length(RMS),o5,'--k');
        title('Radial deviation')
        ax6 = subplot(6,1,6);
        rectangle('Position',[i(6) -0.23 i(2)-1 1.25],'FaceColor',[.95,.95,.95],'Linestyle','none');
        rectangle('Position',[i(7) -0.25 i(2)-1 1.25],'EdgeColor',[.4,.4,.4],'Linewidth',1.2); hold on;
        plot(1:length(RMS),y6_testGPR,'k');
        plot(1:length(RMS),y6_LR,'b');
        plot(1:length(RMS),o6,'--k');
        title('Ulnar deviation')
        xlabel('samples [#]')
        [~, hobj, ~, ~] = legend('GPR','LR','CI(GPR)','Location','Best');
        set(hobj,'linewidth',1.5);
        gca_handles = [ax1,ax2,ax3,ax4,ax5,ax6];
        set(gca_handles,'fontsize',10,'YLim',[-0.25 0.5],'XLim',[0 length(y1_testGPR)+1])
end



% P-values (std) of GPR
% subplot = @(m,n,p) subtightplot(m, n, p, [0.06 0.03], [0.08 0.05], [0.08 0.03]);
% GPR_p=figure('DefaultAxesPosition', [0.1, 0.1, 0.8, 0.8], 'units', 'normalized','position', [0 0 0.39774756441 1],'Color','w');
% ax1=subplot(4,1,1);
% rectangle('Position',[i(3) 0.008 i(2)-1 0.39],'FaceColor',[.95,.95,.95],'Linestyle','none');
% rectangle('Position',[0 0 i(2)-1 0.1],'EdgeColor',[.4,.4,.4],'Linewidth',1.2); hold on;
% plot(1:length(RMS),p1);
% title('Flexion')
% hold on;
% ax2=subplot(4,1,2);
% rectangle('Position',[i(3) 0.008 i(2)-1 0.39],'FaceColor',[.95,.95,.95],'Linestyle','none');
% rectangle('Position',[i(2) 0 i(2)-1 0.1],'EdgeColor',[.4,.4,.4],'Linewidth',1.2); hold on;
% plot(1:length(RMS),p2);
% title('Extension')
% ax3=subplot(4,1,3);
% rectangle('Position',[i(3) 0.008 i(2)-1 0.39],'FaceColor',[.95,.95,.95],'Linestyle','none');
% rectangle('Position',[i(4) 0 i(2)-1 0.1],'EdgeColor',[.4,.4,.4],'Linewidth',1.2); hold on;
% plot(1:length(RMS),p3);
% title('Abduction')
% ylabel('uncertainty [ratio]')
% ax4=subplot(4,1,4);
% rectangle('Position',[i(3) 0.008 i(2)-1 0.39],'FaceColor',[.95,.95,.95],'Linestyle','none');
% rectangle('Position',[i(5) 0 i(2)-3 0.1],'EdgeColor',[.4,.4,.4],'Linewidth',1.2); hold on;
% plot(1:length(RMS),p4);
% title('Adduction')
% xlabel('samples [#]')
% [~, hobj, ~, ~] = legend('p_{GPR}','Location','Best');
% set(hobj,'linewidth',1.5);
% gca_handles = [ax1,ax2,ax3,ax4];
% set(gca_handles,'fontsize',10,'YLim',[0 0.1],'XLim',[0 length(y1_testGPR)+1])

%% Save Figures

clear subplot