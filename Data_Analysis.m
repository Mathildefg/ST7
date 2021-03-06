%clc; clear all; close all;
%% Data analysis - by Tom Baumeister
% The script is written to analyse the obtained data from the fitts law
% test of the Myoband.

subject = 10; % Choose subject number. Each subject has to be done manually to fill the performance matrix.
% Open the directory of the data from the subject
% Load GPR data
load dataGPR.mat
GPR(1).cursor = cursor;                     %set the GPR(1).cursor to have the data from cursor
GPR(1).dof = dof;                           %set the GPR(1).dof to have the data from dof
GPR(1).p = p;                               %set the GPR(1).p to have the data from p
GPR(1).sysOut = sysOut;                     %set the GPR(1).sysOut to have the data from sysOut
GPR(1).target = target;                     %set the GPR(1).target to have the data from target
GPR(1).time = time;                         %set the GPR(1).time to have the data from time
clear cursor dof p sysOut target time       %clear the cursor, dof, p, sysOut, target, time

% Load LR data
load dataLR.mat
LR(1).cursor = cursor;                      %set the LR(1).cursor to have the data from cursor
LR(1).dof = dof;                            %set the LR(1).dof to have the data from dof
LR(1).sysOut = sysOut;                      %set the LR(1).sysOut to have the data from sysOut
LR(1).target = target;                      %set the LR(1).target to have the data from target
LR(1).time = time;                          %set the LR(1).time to have the data from time
clear cursor dof sysOut target time p       %clear the cursor, dof, sysOut, target, time, p

%% Performance measure 1: Time

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

% 47 indexes are expected and are checked. VI SKAL TJEKKE HVOR 
if length(i_lr) == 47               %if the expected length of the array i_lr is 47          
    disp('Yes')                     %display yes
else                                %if the expected length of the array i_lr, is NOT 47, 
    disp('No, there is a mistake')  %display mistake message
end
i_lr = i_lr+1; i_lr = [i_lr;length(LR.time)]; %Add final index

% Show performance: Time it took to reach each target
P = zeros(1,48);                    %create an array of all zeros with dim 1 by 48 called P
for jj = 1:length(P)                %for loop with jj going from 1 through the length of P  %%why is this?? 
    if jj == 1                      %when jj equals 1               %%explain this loop please
        P(jj) = LR.time(i_lr(jj)-1)-1;  
    elseif jj == 48                 
        P(jj) = LR.time(i_lr(48))-LR.time(i_lr(jj-1))-1;  
    else
        P(jj) = LR.time(i_lr(jj)-1)-LR.time(i_lr(jj-1))-1;
    end
end
Score_LR = sum(P);          %Summation of all times is the performance score. 
                            %Smaller number = better (faster time at reaching the target)

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

if length(i_gpr) == 47
    disp('Yes')
else 
    disp('No, there is a mistake')
end
i_gpr = i_gpr+1; i_gpr  = [i_gpr;length(GPR.time)];

% Show performance: Time it took to reach each target
P = zeros(1,48);
for jj = 1:length(P)
    if jj == 1
        P(jj) = GPR.time(i_gpr(jj)-1)-1;
    elseif jj == 48
        P(jj) = GPR.time(i_gpr(48))-GPR.time(i_gpr(jj-1))-1;
    else
        P(jj) = GPR.time(i_gpr(jj)-1)-GPR.time(i_gpr(jj-1))-1;
    end
end
Score_GPR = sum(P);     %Summation of all times is the performance score. 
                        %Smaller number = better (faster time at reaching the target)

% Percentage of GPR being relatively faster. When negative, LR was faster.
B = (Score_LR - Score_GPR)/Score_LR;       

%% Performance measure 2: Path efficiency
% Optimal path and path taken are derived. When cursor was in target, but
% went out of the target before the end of dwell time, the path is included. When
% the cursor went into the target and achieved to stay there for 1 second
% the path taken within the target area is not included. The data points of
% the path were taken at 0.1 sample time. To derive a more realistic
% distance, between each location point the data is interpolated by 100
% data points by use of splines after wich the eucledian distance between
% each new data points is summed.
dt = 0.01;
GPR(1).path.opt = zeros(length(i_gpr),1);       %create an array of zeros with dim length(i_gpr) by 1 called GPR(1).path.opt
LR(1).path.opt = zeros(length(i_lr),1);         %create an array of zeros with dim length(i_lr) by 1 called LR(1).path.opt
GPR(1).path.taken = zeros(length(i_gpr),1);     %create an array of zeros with dim length(i_gpr) by 1 called GPR(1).path.taken
LR(1).path.taken = zeros(length(i_lr),1);       %create an array of zeros with dim length(i_lrr) by 1 called LRR(1).path.taken
for jj = 1:48   %run from jj equal to 1 through 48.    Look into the equations below for Path optimising
    GPR.path.opt(jj) = sqrt(GPR.target(i_gpr(jj)-1,1)^2+GPR.target(i_gpr(jj)-1,2)^2)-0.075; %0.075 is the radius of the target
    LR.path.opt(jj) = sqrt(LR.target(i_lr(jj)-1,1)^2+LR.target(i_lr(jj)-1,2)^2)-0.075;
    if jj == 1  %if jj equals to 1, do:
        gprOut = [0,0;GPR.sysOut(1:i_gpr(jj)-1,1),GPR.sysOut(1:i_gpr(jj)-1,2)];     %GPR out is set to be concatenate of 0,0 and GPR.sysOut(1:i_gpr(jj)-1,1),GPR.sysOut(1:i_gpr(jj)-1,2)
        lrOut = [0,0;LR.sysOut(1:i_lr(jj)-1,1),LR.sysOut(1:i_lr(jj)-1,2)];          %LR out is set to be concatenate of 0,0 and LR.sysOut(1:i_lr(jj)-1,1),LR.sysOut(1:i_lr(jj)-1,2)
        tmax_gpr = length(gprOut);          %returns the length of the longest array dimension in gprOut
        tmax_lr = length(lrOut);            %returns the length of the longest array dimension in lrOut
        t_gpr = [1:dt:tmax_gpr];            %t_gpr = [start:step:stop] -> start at 1, go 0.01 each step, end at tmax_gpr
        t_lr = [1:dt:tmax_lr];              %t_lr = [start:step:stop] -> start at 1, go 0.01 each step, end at tmax_lr
        gprOut_new = interp1(gprOut, t_gpr, 'pchip');   %set gprOut_new to be 1D interpolation of gprOut and t_gpr with 
                                                        %method pchip (Piecewise Cubic Hermite Interpolating Polynomial)
                                                        %connecting the dots
        lrOut_new = interp1(lrOut, t_lr, 'pchip');      %set lrOut_new to be 1D interpolation of lrOut and t_lr with 
                                                        %method pchip (Piecewise Cubic Hermite Interpolating Polynomial)
        count = 2;
        for kk = 3:length(gprOut_new)           %for loop when kk equals 3 through the length of gprOut_new
            if sqrt((GPR.target(i_gpr(jj)-1,1)-gprOut_new(end-(kk-1),1))^2+(GPR.target(i_gpr(jj)-1,2)-gprOut_new(end-(kk-1),2))^2) < 0.075 
                %above: if the ?? is below the radius of the target. ?? dont understand
                count = count + 1;              %add 1 to count, put this as the new count
            else
                break
            end
        end
        real_end_gpr = length(gprOut_new)-count+1;  %real_end_gpr equals the length of gprOut_new - (count+1)
        d_gpr = [gprOut_new(1,1),gprOut_new(1,2);diff(gprOut_new(1:real_end_gpr,:))]; %?? dont understand
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
    elseif jj == 48                 %else if jj is equal to 48 do:
        gprOut = [GPR.sysOut(i_gpr(jj-1):i_gpr(48),1),GPR.sysOut(i_gpr(jj-1):i_gpr(48),2)];
        lrOut = [LR.sysOut(i_lr(jj-1):i_lr(48),1),LR.sysOut(i_lr(jj-1):i_lr(48),2)];
        tmax_gpr = length(gprOut);  %set tmax_gpr to be the length of gprOut
        tmax_lr = length(lrOut);    %set tmax_lr to be the length of lrOut
        t_gpr = [1:dt:tmax_gpr];            %t_gpr = [start:step:stop] -> start at 1, go 0.01 each step, end at tmax_gpr
        t_lr = [1:dt:tmax_lr];              %t_lr = [start:step:stop] -> start at 1, go 0.01 each step, end at tmax_lr
        gprOut_new = interp1(gprOut, t_gpr, 'pchip');   %set gprOut_new to be 1D interpolation of gprOut and t_gpr with 
                                                        %method pchip (Piecewise Cubic Hermite Interpolating Polynomial)
                                                        %connecting the dots
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
        lrOut = [LR.sysOut(i_lr(jj-1):i_lr(jj)-1,1),LR.sysOut(i_lr(jj-1):i_lr(jj)-1,2)];
        tmax_gpr = length(gprOut);
        tmax_lr = length(lrOut);
        t_gpr = [1:dt:tmax_gpr];
        t_lr = [1:dt:tmax_lr];
        gprOut_new = interp1(gprOut, t_gpr, 'pchip');
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
    LR.path.taken(jj) = sum(sqrt(sum(d_lr.*d_lr,2)));
end
clear count d_gpr d_lr dt gprOut gprOut_new jj kk lrOut lrOut_new real_end_gpr real_end_lr t_gpr t_lr tmax_gpr tmax_lr
% Path efficiency calculated
path_ratio_gpr=(GPR.path.taken - GPR.path.opt)./GPR.path.opt;
path_ratio_lr=(LR.path.taken - LR.path.opt)./LR.path.opt;

%% Diagonal and non-diagonal path efficiencies
% Find the indexes of diagonal and non-diagonal
% GPR
gpr_targets = GPR.target([1;i_gpr(1:end-1)],:);
gpr_targets_check = zeros(length(gpr_targets),1);
for jj = 1:length(gpr_targets)
    if abs(gpr_targets(jj,1)) == abs(gpr_targets(jj,2))
        gpr_targets_check(jj) = 1;
    end
end
i_gpr_d = find(gpr_targets_check);          %find the indexes and values of the nonzero elements in gpr_targets_check
i_gpr_nd = find(gpr_targets_check==0);      %find the index equal to 0 in the gpr_targets_check

% LR
lr_targets = LR.target([1;i_lr(1:end-1)],:);
lr_targets_check = zeros(length(lr_targets),1);
for jj = 1:length(lr_targets)
    if abs(lr_targets(jj,1)) == abs(lr_targets(jj,2))
        lr_targets_check(jj) = 1;
    end
end
i_lr_d = find(lr_targets_check);            %find the indexes and values of the nonzero elements in lr_targets_check
i_lr_nd = find(lr_targets_check==0);        %find the index equal to 0 in the lr_targets_check

% Path ratios of diagonal and non-diagonal targets for lr and gpr.
path_ratio_gpr_d=(GPR.path.taken(i_gpr_d) - GPR.path.opt(i_gpr_d))./GPR.path.opt(i_gpr_d);
path_ratio_gpr_nd=(GPR.path.taken(i_gpr_nd) - GPR.path.opt(i_gpr_nd))./GPR.path.opt(i_gpr_nd);
path_ratio_lr_d=(LR.path.taken(i_lr_d) - LR.path.opt(i_lr_d))./LR.path.opt(i_lr_d);
path_ratio_lr_nd=(LR.path.taken(i_lr_nd) - LR.path.opt(i_lr_nd))./LR.path.opt(i_lr_nd);

%% Performance matrices: total and diagonal vs non-diagonal
% All performance measures are put into 1 matrix. This is done manually for
% each subject.

%Performance = zeros(10,6);
Performance(subject,1) = Score_LR;
Performance(subject,2) = Score_GPR;
Performance(subject,3) = B;
Performance(subject,4) = mean(path_ratio_lr);
Performance(subject,5) = mean(path_ratio_gpr);
Performance(subject,6) = (Performance(subject,4)-Performance(subject,5))./Performance(subject,4);

%Performance_d = zeros(10,6);
Performance_d(subject,1) = mean(path_ratio_gpr_d);
Performance_d(subject,2) = mean(path_ratio_gpr_nd);
Performance_d(subject,3) = mean(path_ratio_lr_d);
Performance_d(subject,4) = mean(path_ratio_lr_nd);
Performance_d(subject,5) = (Performance_d(subject,3)-Performance(subject,1))./Performance(subject,3);
Performance_d(subject,6) = (Performance_d(subject,4)-Performance(subject,2))./Performance(subject,4);

%% EXTRAS: This include trajectory simultion and a sample time check
% % Simulate Trajectory 
% close all
% hfig = figure;
% set(hfig, 'units', 'normalized','position', [0 0 1 1],'menubar', 'none','renderer','painters','name','Experiment','numbertitle','off')
% % Add target marker
% htarget = plot(LR.target(i_lr(jj)-1,1), LR.target(i_lr(jj)-1,2), 'ro','markersize',50,'linewidth', 2, 'buttondownfcn', 'uiresume'); hold on
% % Add position (of controlled system) marker
% hpos = plot(0, 0, 'ko', 'markersize', 12, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0.1 0.1 0.1], 'linewidth', 2, 'buttondownfcn', 'uiresume');
% axis([-1 1 -1 1])
% set(gca, 'color', [.98,.98,.98],'DataAspectRatio',[1 1 1], 'units', 'normalized', 'position', [0 0 1 1], 'xtick', 0, 'ytick', 0, 'LineWidth', 3,'GridColor','k')%,'drawmode', 'fast')
% grid on
% % Wait for start. When clicked on the cursor the simulation starts
% uiwait
% % Hide mouse pointer
% set(gcf, 'pointer', 'custom','pointershapecdata', nan(16,16))
% drawnow
% 
% x1 = LR.sysOut(i_lr(jj-1):i_lr(jj)-1,1);
% y1 = LR.sysOut(i_lr(jj-1):i_lr(jj)-1,2);
% % plot(x1,y1,'bx');
% 
% dt = 0.1;
% tic
% for ii = 1:length(x1)
%     set(hpos, 'xdata', x1(ii), 'ydata', y1(ii)); hold on;
%     plot(x1(ii),y1(ii),'bx')
%     drawnow  
%     while (toc < (ii) * dt)
%     end
% end
% 
% %% Sample Time check  (cut in pieces to make it fit on page instead of
% all the way to the right 
%FIRST PART:
% dt_gpr = diff(GPR.time(i(1):i(2)));diff(GPR.time(i(2)+1:i(3)));
%diff(GPR.time(i(3)+1:i(4)));diff(GPR.time(i(4)+1:i(5)));
%diff(GPR.time(i(5)+1:i(6)));diff(GPR.time(i(6)+1:i(7)));
%diff(GPR.time(i(7)+1:i(8)));diff(GPR.time(i(8)+1:i(9)));
%diff(GPR.time(i(9)+1:i(10)));diff(GPR.time(i(10)+1:i(11)));
%diff(GPR.time(i(11)+1:i(12)));diff(GPR.time(i(12)+1:i(13)));
%diff(GPR.time(i(13)+1:i(14)));diff(GPR.time(i(14)+1:i(15)));
%diff(GPR.time(i(15)+1:i(16)));diff(GPR.time(i(16)+1:i(17)));
%diff(GPR.time(i(17)+1:i(18)));diff(GPR.time(i(18)+1:i(19)));
%diff(GPR.time(i(19)+1:i(20)));
%SECOND PART:
% dt_lr = diff(LR.time(i(1):i(2)));diff(LR.time(i(2)+1:i(3)));
%diff(LR.time(i(3)+1:i(4)));diff(LR.time(i(4)+1:i(5)));
%diff(LR.time(i(5)+1:i(6)));diff(LR.time(i(6)+1:i(7)));
%diff(LR.time(i(7)+1:i(8)));diff(LR.time(i(8)+1:i(9)));
%diff(LR.time(i(9)+1:i(10)));diff(LR.time(i(10)+1:i(11)));
%diff(LR.time(i(11)+1:i(12)));diff(LR.time(i(12)+1:i(13)));
%diff(LR.time(i(13)+1:i(14)));diff(LR.time(i(14)+1:i(15)));
%diff(LR.time(i(15)+1:i(16)));diff(LR.time(i(16)+1:i(17)));
%diff(LR.time(i(17)+1:i(18)));diff(LR.time(i(18)+1:i(19)));diff(LR.time(i(19)+1:i(20)));
% 
% disp('The mean dt for GPR and LR are: ',mean(dt_gpr),' ',mean(dt_lr));
% disp('The std of the dt for GPR and LR are: ',std(dt_gpr),' ',std(dt_lr));

