%% PRØVE

fileFolder = dir(uigetdir);
fileFolder = fileFolder.folder;
addpath(fileFolder);

% Get a list of all .mat-files in the folder.
listOfFiles = dir(fullfile(fileFolder, '*.mat'));
% Load the files
EMG= []; MVC=[]; RMS= []; 

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
            
            
        end