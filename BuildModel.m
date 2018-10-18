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
EMG= []; RMS= []; GeneratedProfile = [];movement=[];
for k = 1 : length(listOfFiles)

  load(listOfFiles(k).name, 'emg');
  data{k} = emg.values;
  EMG = [EMG; data{1,k}];
  
  RMS_temp = [];
  for i=1 : 8
      RMS_temp(:,i) = transpose(my_rms(data{1,k}(:,i),40,30,0));
     
  end
  RMS = [RMS; RMS_temp];
  
  GeneratedProfile = [GeneratedProfile my_rms(data{1,k}, 40, 30, 0)];
  
  %To be used for indexing - the lengths of the RMS from each movement:
  f(k) = length(RMS_temp);
  
end
clear emg RMS_temp;

%Assigning a number for each movement to be used in indexing. Feasible
%because the files are always loaded alphabetically9
movement(1:sum(f(1:3))) =2;
movement(sum(f(1:3))+1: sum(f(1:6)))=1;
movement(sum(f(1:6))+1:sum(f(1:9)))=6;
movement(sum(f(1:9))+1:sum(f(1:12)))=3;
movement(sum(f(1:12))+1:sum(f(1:15)))=7;
movement(sum(f(1:15))+1:sum(f(1:19)))= 5;
movement(sum(f(1:19))+1:sum(f(1:21)))=4;

%The index vector
[dofs, i, ~] = unique(movement);
i = sort(i);


figure(1)
subplot(3,1,1)
plot(EMG)
subplot(3,1,2)
plot(RMS)
subplot(3,1,3)
plot(GeneratedProfile)


%% Arrangement in predictors and target

%Predictor - RMS of each EMG-channel of the myo armband for 25%, 50%, 75% of MVC. 
x_flex = RMS(i(1):i(2)-1,:);  %Vi trækker 1 fra fordi i(2) er første sample i næste bevægelse.
x_ext = RMS(i(2):i(3)-1,:);
x_rd = RMS(i(3):i(4)-1,:);
x_ud = RMS(i(4):i(5)-1,:);
x_pro = RMS(i(5):i(6)-1,:);
x_sup = RMS(i(6):end,:);

%Target values - The generated profile. 
% y_flex
% y_ext
% y_rd
% y_ud
% y_pro
% y_sup