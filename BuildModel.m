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
EMG= []; RMS= []; GeneratedProfile = [];
for k = 1 : length(listOfFiles)
  %data{k}= eval(['load ' listOfFiles(k).name]);
  %sprintf('data%d', k) = load(listOfFiles(k).name);
  %data = 
  %filename = sprintf('data%d', k);
 
  load(listOfFiles(k).name, 'emg');
  data{k} = emg.values;
  EMG = [EMG; data{1,k}];
  
  RMS_temp = [];
  for i=1 : 8
      RMS_temp(:,i) = transpose(my_rms(data{1,k}(:,i),20,15,0));
  end
  RMS = [RMS; RMS_temp];
 
  
  GeneratedProfile = [GeneratedProfile my_rms(data{1,k}, 20, 15, 0)];
end
clear emg;

figure(1)
subplot(3,1,1)
plot(EMG)
subplot(3,1,2)
plot(RMS)
subplot(3,1,3)
plot(GeneratedProfile)


%% Feature extraction and arrangement in predictors and target

%Predictor - RMS of each EMG-channel of the myo armband for 25%, 50%, 75% of MVC. 
% x_flex
% x_ext
% x_rd 
% x_ud
% x_pro
% x_sup
% 
% %Target values - The generated profile. 
% y_flex
% y_ext
% y_rd
% y_ud
% y_pro
% y_sup