%% IMPORT DATA
for Testperson=1:10
%Testperson= 1
fileFolder = dir(uigetdir);
fileFolder = fileFolder.folder;
addpath(fileFolder);

%GPR
load('dataGPR_mc.mat');
GPR.cursor = cursor;                     
GPR.dof = dof;                           
GPR.p = p;                              
GPR.sysOut = sysOut;                     
GPR.target = target;                     
GPR.time = time; 
GPR.ci = ci; 

%GPRc
load('dataGPR_uc.mat');
GPRc.cursor = cursor;                     
GPRc.dof = dof;                           
GPRc.p = p;                              
GPRc.sysOut = sysOut;                     
GPRc.target = target;                     
GPRc.time = time; 
GPRc.ci = ci; 

%LR
load('dataLR.mat');
LR.cursor = cursor;                     
LR.dof = dof;                           
LR.p = p;                              
LR.sysOut = sysOut;                     
LR.target = target;                     
LR.time = time; 

%% INDEX of difficulty

fileID = fopen('Fitt_real_2dof.txt','r');           
Target = fscanf(fileID,'%f %f %f %f',[4 inf])';



ID = [];
W= 25; %Target-width in pixels 
% The -1 to 1 coordinatesystem is normalized by the following:
pix_width = 1368; 
pix_height = 891.5;

for i = 1 : length(Target)
    %For at få target-koordinaterne i pixels denormaliseres de:
x_pix= Target(i,1)*pix_width;
y_pix = Target(i,2)*pix_height;
d(i)= sqrt(x_pix^2+y_pix^2);
ID(i) = log2(d(i)/W + 1);
end

% For d_unique
% d_unique = unique(d);
% 
% for k = 1 : length(d_unique)
%     %For at få target-koordinaterne i pixels denormaliseres de:
% 
% ID_unique(k) = log2(d_unique(k)/W + 1);
% end

%% MOVEMENT TIME
% LR
% Get the indexes for each trial
i_lr = zeros(100,1);                %i_LR. create an array of all zeros of dimension 100 by 1.
j = 1;                              %j = one
for ii = 1:length(LR.time)-1        %from one through length of the LR.time - 1.  %%Why this iteration??
    if LR.time(ii+1)-LR.time(ii)>1  %After each trial there is a pause of 1 second that is used to find the indexes
        A = LR.time(ii);            %
        i_lr(j) = find(LR.time==A); %
        j=j+1;                      %
    end    
end
i_lr = nonzeros(i_lr);              %get the nonzero elements in i_lr, put them as i_lr

% 31 indexes are expected and are checked. 
if length(i_lr) == 31               %if the expected length of the array i_lr is 47          
    disp('Yes')                     %display yes
else                                %if the expected length of the array i_lr, is NOT 47, 
    disp('No, there is a mistake')  %display mistake message
end
i_lr = i_lr+1; i_lr = [i_lr;length(LR.time)]; %Add final index

% Show time it took to reach each target
trial_time_LR = zeros(1,32);                    %create an array of all zeros with dim 1 by 32 called P
for jj = 1:length(trial_time_LR)                %for loop with jj going from 1 through the length of P 
    if jj == 1                      %when jj equals 1         
        trial_time_LR(jj) = LR.time(i_lr(jj)-1)-1;  
    elseif jj == 32                 
        trial_time_LR(jj) = LR.time(i_lr(32))-LR.time(i_lr(jj-1))-1;  
    else
        trial_time_LR(jj) = LR.time(i_lr(jj)-1)-LR.time(i_lr(jj-1))-1;
    end
end

fails_LR=0;
for g=1:length(trial_time_LR)
    if trial_time_LR(g)>=15  
       trial_time_LR(g)=15;
       fails_LR = fails_LR+1;
    end
end

TotalTime_LR = sum(trial_time_LR);

% GPR (Same is done as in LR)
i_gpr = zeros(100,1);
j = 1;
for ii = 1:length(GPR.time)-1
    if GPR.time(ii+1)-GPR.time(ii)>1
        A = GPR.time(ii);
        i_gpr(j) = find(GPR.time==A);
        j=j+1;
    end
end
i_gpr = nonzeros(i_gpr);

if length(i_gpr) == 31
    disp('Yes')
else 
    disp('GPR: No, there is a mistake')
end
i_gpr = i_gpr+1; i_gpr  = [i_gpr;length(GPR.time)];

% Show performance: Time it took to reach each target
trial_time_GPR = zeros(1,32);
for jj = 1:length(trial_time_GPR)
    if jj == 1
        trial_time_GPR(jj) = GPR.time(i_gpr(jj)-1)-1;
    elseif jj == 32
        trial_time_GPR(jj) = GPR.time(i_gpr(32))-GPR.time(i_gpr(jj-1))-1;
    else
        trial_time_GPR(jj) = GPR.time(i_gpr(jj)-1)-GPR.time(i_gpr(jj-1))-1;
    end
end

fails_GPR=0;
for g=1:length(trial_time_GPR)
    if trial_time_GPR(g)>=15  
       trial_time_GPR(g)=15;
       fails_GPR = fails_GPR+1;
    end
end

TotalTime_GPR = sum(trial_time_GPR);

%GPRc
i_gpr_c = zeros(100,1);
j = 1;
for ii = 1:length(GPRc.time)-1
    if GPRc.time(ii+1)-GPRc.time(ii)>1
        A = GPRc.time(ii);
        i_gpr_c(j) = find(GPRc.time==A);
        j=j+1;
    end
end
i_gpr_c = nonzeros(i_gpr_c);

if length(i_gpr_c) == 31
    disp('Yes')
else 
    disp('GPRc: No, there is a mistake')
end
i_gpr_c = i_gpr_c+1; i_gpr_c  = [i_gpr_c;length(GPRc.time)];

% Show performance: Time it took to reach each target
trial_time_GPRc = zeros(1,32);
for jj = 1:length(trial_time_GPRc)
    if jj == 1
        trial_time_GPRc(jj) = GPRc.time(i_gpr_c(jj)-1)-1;
    elseif jj == 32
        trial_time_GPRc(jj) = GPRc.time(i_gpr_c(32))-GPRc.time(i_gpr_c(jj-1))-1;
    else
        trial_time_GPRc(jj) = GPRc.time(i_gpr_c(jj)-1)-GPRc.time(i_gpr_c(jj-1))-1;
    end
end

fails_GPRc=0;
for g=1:length(trial_time_GPRc)
    if trial_time_GPRc(g)>=15  
       trial_time_GPRc(g)=15;
       fails_GPRc = fails_GPRc+1;
    end
end

TotalTime_GPRc = sum(trial_time_GPRc); 
%% Calculate THROUGHPUT
%Throughput = zeros(10, 3);
Throughput = [mean(ID/trial_time_LR),mean(ID/trial_time_GPR),mean(ID/trial_time_GPRc)]; 


%% Path Efficiency
% Optimal path and path taken are derived. When cursor was in target, but
% went out of the target before the end of dwell time, the path is included. When
% the cursor went into the target and achieved to stay there for 1 second
% the path taken within the target area is not included. The data points of
% the path were taken at 0.1 sample time. To derive a more realistic
% distance, between each location point the data is interpolated by 100
% data points by use of splines after which the eucledian distance between
% each new data points is summed.
dt = 0.01;
GPR.path.opt = zeros(length(i_gpr),1);       
LR.path.opt = zeros(length(i_lr),1);  
GPRc.path.opt = zeros(length(i_gpr_c),1);


GPR.path.taken = zeros(length(i_gpr),1);     
GPRc.path.taken = zeros(length(i_gpr_c),1);  
LR.path.taken = zeros(length(i_lr),1);      

for jj = 1:32   
    %Optimal path calculated by pythagoras
    GPR.path.opt(jj) = sqrt(GPR.target(i_gpr(jj)-1,1)^2+GPR.target(i_gpr(jj)-1,2)^2)-0.075;  %0.075 is the radius of the target
    GPRc.path.opt(jj) = sqrt(GPRc.target(i_gpr_c(jj)-1,1)^2+GPRc.target(i_gpr_c(jj)-1,2)^2)-0.075;
    LR.path.opt(jj) = sqrt(LR.target(i_lr(jj)-1,1)^2+LR.target(i_lr(jj)-1,2)^2)-0.075;
   
    
    if jj == 1  %if jj equals to 1, do:
        gprOut = [0,0;GPR.sysOut(1:i_gpr(jj)-1,1),GPR.sysOut(1:i_gpr(jj)-1,2)];     %GPR out is set to be concatenate of 0,0 and GPR.sysOut(1:i_gpr(jj)-1,1),GPR.sysOut(1:i_gpr(jj)-1,2)
        gprOut_c = [0,0;GPRc.sysOut(1:i_gpr_c(jj)-1,1),GPRc.sysOut(1:i_gpr_c(jj)-1,2)];
        lrOut = [0,0;LR.sysOut(1:i_lr(jj)-1,1),LR.sysOut(1:i_lr(jj)-1,2)];          %LR out is set to be concatenate of 0,0 and LR.sysOut(1:i_lr(jj)-1,1),LR.sysOut(1:i_lr(jj)-1,2)
        
        tmax_gpr = length(gprOut);          %returns the length of the longest array dimension in gprOut
        tmax_gpr_c = length(gprOut_c); 
        tmax_lr = length(lrOut);            %returns the length of the longest array dimension in lrOut
        
        t_gpr = [1:dt:tmax_gpr];            %t_gpr = [start:step:stop] -> start at 1, go 0.01 each step, end at tmax_gpr
        t_gpr_c = [1:dt:tmax_gpr_c];
        t_lr = [1:dt:tmax_lr];              %t_lr = [start:step:stop] -> start at 1, go 0.01 each step, end at tmax_lr
        
        gprOut_new = interp1(gprOut, t_gpr, 'pchip');   %set gprOut_new to be 1D interpolation of gprOut and t_gpr with 
                                                        %method pchip (Piecewise Cubic Hermite Interpolating Polynomial)
        gprOut_c_new = interp1(gprOut_c, t_gpr_c, 'pchip');                           %connecting the dots
        lrOut_new = interp1(lrOut, t_lr, 'pchip');      %set lrOut_new to be 1D interpolation of lrOut and t_lr with 
                                                        %method pchip (Piecewise Cubic Hermite Interpolating Polynomial)
        count = 2;
        for kk = 3:length(gprOut_new)          
            if sqrt((GPR.target(i_gpr(jj)-1,1)-gprOut_new(end-(kk-1),1))^2+(GPR.target(i_gpr(jj)-1,2)-gprOut_new(end-(kk-1),2))^2) < 0.075 
             
                count = count + 1;              %add 1 to count, put this as the new count
            else
                break
            end
        end
        real_end_gpr = length(gprOut_new)-count+1;  %real_end_gpr equals the length of gprOut_new - (count+1)
        d_gpr = [gprOut_new(1,1),gprOut_new(1,2);diff(gprOut_new(1:real_end_gpr,:))]; 
        
        count = 2;
        for kk = 3:length(gprOut_c_new)          
            if sqrt((GPRc.target(i_gpr_c(jj)-1,1)-gprOut_c_new(end-(kk-1),1))^2+(GPRc.target(i_gpr_c(jj)-1,2)-gprOut_c_new(end-(kk-1),2))^2) < 0.075 
                %above: if the ?? is below the radius of the target. ?? dont understand
                count = count + 1;              %add 1 to count, put this as the new count
            else
                break
            end
        end
        real_end_gpr_c = length(gprOut_c_new)-count+1;  %real_end_gpr equals the length of gprOut_new - (count+1)
        d_gpr_c = [gprOut_c_new(1,1),gprOut_c_new(1,2);diff(gprOut_c_new(1:real_end_gpr_c,:))]; %?? dont understand
        
        
        count = 2;                              %set count to be equal to 2
        for kk = 3:length(lrOut_new)            %for loop when kk equals 3 through length of lrOut_new
            if sqrt((LR.target(i_lr(jj)-1,1)-lrOut_new(end-(kk-1),1))^2+(LR.target(i_lr(jj)-1,2)-lrOut_new(end-(kk-1),2))^2) < 0.075
                %above: if the ?? is below the radius of the target. ?? dont understand
                count = count + 1;              %add 1 to count, put this as the new count
            else
                break
            end
        end
        real_end_lr = length(lrOut_new)-count+1;    %real_end_lr equals the length of lrOut_new - (count+1)
        d_lr = [lrOut_new(1,1),lrOut_new(1,2);diff(lrOut_new(1:real_end_lr,:))]; %?? dont understand
        
    elseif jj == 32                 %else if jj is equal to 48 do:
        gprOut = [GPR.sysOut(i_gpr(jj-1):i_gpr(32),1),GPR.sysOut(i_gpr(jj-1):i_gpr(32),2)];
        gprOut_c = [GPRc.sysOut(i_gpr_c(jj-1):i_gpr_c(32),1),GPRc.sysOut(i_gpr_c(jj-1):i_gpr_c(32),2)];
        lrOut = [LR.sysOut(i_lr(jj-1):i_lr(32),1),LR.sysOut(i_lr(jj-1):i_lr(32),2)];
        
        tmax_gpr = length(gprOut);  %set tmax_gpr to be the length of gprOut
         tmax_gpr_c = length(gprOut_c);
        tmax_lr = length(lrOut);    %set tmax_lr to be the length of lrOut
        
        t_gpr = [1:dt:tmax_gpr];            %t_gpr = [start:step:stop] -> start at 1, go 0.01 each step, end at tmax_gpr
        t_gpr_c = [1:dt:tmax_gpr_c];
        t_lr = [1:dt:tmax_lr];              %t_lr = [start:step:stop] -> start at 1, go 0.01 each step, end at tmax_lr
        gprOut_new = interp1(gprOut, t_gpr, 'pchip');   %set gprOut_new to be 1D interpolation of gprOut and t_gpr with 
                                                        %method pchip (Piecewise Cubic Hermite Interpolating Polynomial)
        gprOut_c_new = interp1(gprOut_c, t_gpr_c, 'pchip');                                                %connecting the dots
        lrOut_new = interp1(lrOut, t_lr, 'pchip');      %set lrOut_new to be 1D interpolation of lrOut and t_lr with 
                                                        %method pchip (Piecewise Cubic Hermite Interpolating Polynomial)
        count = 2;
        for kk = 3:length(gprOut_new)
            if sqrt((GPR.target(i_gpr(jj)-1,1)-gprOut_new(end-(kk-1),1))^2+(GPR.target(i_gpr(jj)-1,2)-gprOut_new(end-(kk-1),2))^2) < 0.075
                count = count + 1;
            else
                break
            end
        end
        real_end_gpr = length(gprOut_new)-count+1;
        d_gpr = [gprOut_new(1,1),gprOut_new(1,2);diff(gprOut_new(1:real_end_gpr,:))];
        
        count = 2;
        for kk = 3:length(gprOut_c_new)
            if sqrt((GPRc.target(i_gpr_c(jj)-1,1)-gprOut_c_new(end-(kk-1),1))^2+(GPRc.target(i_gpr_c(jj)-1,2)-gprOut_c_new(end-(kk-1),2))^2) < 0.075
                count = count + 1;
            else
                break
            end
        end
        real_end_gpr_c = length(gprOut_c_new)-count+1;
        d_gpr_c = [gprOut_c_new(1,1),gprOut_c_new(1,2);diff(gprOut_c_new(1:real_end_gpr_c,:))];
        
        count = 2;
        for kk = 3:length(lrOut_new)
            if sqrt((LR.target(i_lr(jj)-1,1)-lrOut_new(end-(kk-1),1))^2+(LR.target(i_lr(jj)-1,2)-lrOut_new(end-(kk-1),2))^2) < 0.075
                count = count + 1;
            else
                break
            end
        end
        real_end_lr = length(lrOut_new)-count+1;
        d_lr = [lrOut_new(1,1),lrOut_new(1,2);diff(lrOut_new(1:real_end_lr,:))];
    else
        gprOut = [GPR.sysOut(i_gpr(jj-1):i_gpr(jj)-1,1),GPR.sysOut(i_gpr(jj-1):i_gpr(jj)-1,2)];
        gprOut_c = [GPRc.sysOut(i_gpr_c(jj-1):i_gpr_c(jj)-1,1),GPRc.sysOut(i_gpr_c(jj-1):i_gpr_c(jj)-1,2)];
        lrOut = [LR.sysOut(i_lr(jj-1):i_lr(jj)-1,1),LR.sysOut(i_lr(jj-1):i_lr(jj)-1,2)];
       
        tmax_gpr = length(gprOut);
        tmax_gpr_c = length(gprOut_c);
        tmax_lr = length(lrOut);
        t_gpr = [1:dt:tmax_gpr];
        t_gpr_c = [1:dt:tmax_gpr_c];
        t_lr = [1:dt:tmax_lr];
        
        gprOut_new = interp1(gprOut, t_gpr, 'pchip');
        gprOut_c_new = interp1(gprOut_c, t_gpr_c, 'pchip');
        lrOut_new = interp1(lrOut, t_lr, 'pchip');
        
        count = 2;
        for kk = 3:length(gprOut_new)
            if sqrt((GPR.target(i_gpr(jj)-1,1)-gprOut_new(end-(kk-1),1))^2+(GPR.target(i_gpr(jj)-1,2)-gprOut_new(end-(kk-1),2))^2) < 0.075
                count = count + 1;
            else
                break
            end
        end
        real_end_gpr = length(gprOut_new)-count+1;
        d_gpr = [gprOut_new(1,1),gprOut_new(1,2);diff(gprOut_new(1:real_end_gpr,:))];
        
        count = 2;
        for kk = 3:length(gprOut_c_new)
            if sqrt((GPRc.target(i_gpr_c(jj)-1,1)-gprOut_c_new(end-(kk-1),1))^2+(GPRc.target(i_gpr_c(jj)-1,2)-gprOut_c_new(end-(kk-1),2))^2) < 0.075
                count = count + 1;
            else
                break
            end
        end
        real_end_gpr_c = length(gprOut_c_new)-count+1;
        d_gpr_c = [gprOut_c_new(1,1),gprOut_c_new(1,2);diff(gprOut_c_new(1:real_end_gpr_c,:))];
        
        count = 2;
        for kk = 3:length(lrOut_new)
            if sqrt((LR.target(i_lr(jj)-1,1)-lrOut_new(end-(kk-1),1))^2+(LR.target(i_lr(jj)-1,2)-lrOut_new(end-(kk-1),2))^2) < 0.075
                count = count + 1;
            else
                break
            end
        end
        real_end_lr = length(lrOut_new);
        d_lr = [lrOut_new(1,1),lrOut_new(1,2);diff(lrOut_new(1:real_end_lr,:))];
    end
    GPR.path.taken(jj) = sum(sqrt(sum(d_gpr.*d_gpr,2)));
    GPRc.path.taken(jj) = sum(sqrt(sum(d_gpr_c.*d_gpr_c,2)));
    LR.path.taken(jj) = sum(sqrt(sum(d_lr.*d_lr,2)));
   
    if trial_time_GPR(jj)>=15 
        GPR.path.opt(jj)=0;
    end
    
    if trial_time_GPRc(jj)>=15 
        GPRc.path.opt(jj)= 0;
    end
    
    if trial_time_LR(jj)>=15
        LR.path.opt(jj)=0;
    end
    
    
end
clear count d_gpr d_lr dt gprOut gprOut_new jj kk lrOut lrOut_new real_end_gpr real_end_lr t_gpr t_lr tmax_gpr tmax_lr

% Path efficiency calculated
% path_ratio_gpr=(GPR.path.taken - GPR.path.opt)./GPR.path.opt;
% path_ratio_gpr_c=(GPRc.path.taken - GPRc.path.opt)./GPRc.path.opt;
% path_ratio_lr=(LR.path.taken - LR.path.opt)./LR.path.opt;

    
path_ratio_gpr=(GPR.path.opt./GPR.path.taken)*100;
path_ratio_gpr_c=(GPRc.path.opt./GPRc.path.taken)*100;
path_ratio_lr=(LR.path.opt./LR.path.taken)*100;



PathEfficiency = [mean(path_ratio_lr) mean(path_ratio_gpr) mean(path_ratio_gpr_c)];

%% COMPLETION RATE

CompletionRate = [((32-fails_LR)/32)*100 ((32-fails_GPR)/32)*100 ((32-fails_GPRc)/32)*100]

%% SAMLES TIL ÉN MATRIX
%FittsPerformanceMetrics= zeros(10,9);
FittsPerformanceMetrics(Testperson,:) = [Throughput,PathEfficiency,CompletionRate]
end