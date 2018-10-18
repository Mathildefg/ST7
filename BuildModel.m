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
EMG= []; RMS = [];
for k = 1 : length(listOfFiles)
  %data{k}= eval(['load ' listOfFiles(k).name]);
  %sprintf('data%d', k) = load(listOfFiles(k).name);
  %data = 
  %filename = sprintf('data%d', k);
 
  load(listOfFiles(k).name, 'emg');
  data{k} = emg.values;
  EMG = [EMG; data{1,k}];
  RMS = [RMS rms(data{1,k}, 20, 15, 0)];
end
clear emg;

% GeneratedProfile = 

figure(1)
title('EMG-data, RMS og Generated profile');
subplot(2,1,1)
plot(EMG)
subplot(2,1,2)
plot(RMS)


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