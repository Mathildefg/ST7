clc; clear; close all;
%% ST7 - Gaussian Process Regression (GPR)
%This script uses data collected using the MATLAB-app
%'CollectTrainingData', to build a GPR-model.

%% Loading and structuring data
% Prompt the user to choose the folder containing the data of the
% respective test subject. 
fileFolder = dir(uigetdir);
fileFolder = fileFolder.folder;

% Check to make sure that folder actually exists.  Warn user if it doesn't.
if ~isdir(fileFolder.folder)
  errorMessage = sprintf('Error: The following folder does not exist:\n%s', fileFolder);
  uiwait(warndlg(errorMessage));
  return;
end

% Get a list of all .mat-files in the folder.
listOfFiles = dir(fullfile(fileFolder, '*.mat'));
% Load the files 
for k = 1 : length(listOfFiles)
  eval(['load ' listOfFiles(k).name]);
end



% EMG
% RMS= 
% GeneratedProfile = 




%% Feature extraction and arrangement in predictors and target

%Predictor - RMS of each EMG-channel of the myo armband for 25%, 50%, 75% of MVC. 
x_flex
x_ext
x_rd 
x_ud
x_pro
x_sup

%Target values - The generated profile. 
y_flex
y_ext
y_rd
y_ud
y_pro
y_sup