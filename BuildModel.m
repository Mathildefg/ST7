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
  
  %To be used for indexing:
  f(k) = length(RMS_temp);
  
  %Assigning a number for each movement to be used in indexing
  switch listOfFiles(k).name
      case {'Flexion_1.mat', 'Flexion_2.mat', 'Flexion_3.mat'}
          temp_array(1:f(k))= 1;
          movement = [movement temp_array];
          
      case {'Extension_1.mat', 'Extension_2.mat', 'Extension_3.mat'}
          temp_array(1:f(k))= 2;
          movement = [movement temp_array];
          
      case {'Radial_deviation_1.mat', 'Radial_deviation_2.mat', 'Radial_deviation_3.mat'}
          temp_array(1:f(k))= 3;
          movement = [movement temp_array];
          
      case {'Ulnar_deviation_1.mat', 'Ulnar_deviation_2.mat', 'Ulnar_deviation_3.mat'}
          temp_array(1:f(k))= 4;
          movement = [movement temp_array];
          
      case {'Supination_1.mat', 'Supination_2.mat', 'Supination_3.mat'}
          temp_array(1:f(k))= 5;
          movement = [movement temp_array];
          
      case {'Pronation_1.mat', 'Pronation_2.mat', 'Pronation_3.mat'}
          temp_array(1:f(k))= 6;
          movement = [movement temp_array];
          
      case {'Rest_1.mat', 'Rest_2.mat', 'Rest_3.mat'}
          temp_array(1:f(k))= 7;
          movement = [movement temp_array];
          
      otherwise
          fprintf('Unknown file included\n')
  end
  
end
clear emg RMS_temp;

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