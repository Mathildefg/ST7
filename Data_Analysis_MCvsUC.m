%clc; clear all; close all;
%% Data analysis - by Tom Baumeister
% The script is written to analyse the obtained data from the fitts law
% test of the Myoband.

subject = 12; % Choose subject number. Each subject has to be done manually to fill the performance matrix.
% Open the directory of the data from the subject
% Load GPR_uc data
dataLoad = dir(uigetdir);
dataLoad = dataLoad.folder;
addpath(dataLoad);
load dataGPR_uc.mat
GPR_uc(1).cursor = cursor;                     %set the GPR(1).cursor to have the data from cursor
GPR_uc(1).dof = dof;                           %set the GPR(1).dof to have the data from dof
GPR_uc(1).p = p;                               %set the GPR(1).p to have the data from p
GPR_uc(1).sysOut = sysOut;                     %set the GPR(1).sysOut to have the data from sysOut
GPR_uc(1).target = target;                     %set the GPR(1).target to have the data from target
GPR_uc(1).time = time;                         %set the GPR(1).time to have the data from time
clear cursor dof p sysOut target time       %clear the cursor, dof, p, sysOut, target, time

% Load GPR_mc data
load dataGPR_mc.mat
GPR_mc(1).cursor = cursor;
GPR_mc(1).dof = dof;
GPR_mc(1).p = p;
GPR_mc(1).sysOut = sysOut;
GPR_mc(1).target = target;
GPR_mc(1).time = time;
clear cursor dof sysOut target time p

%% Performance measure 1: Time

% GPR_mc
% Get the indexes for each trial
i_GPR_mc = zeros(100,1);                %i_GPR_mc. create an array of all zeros of dimension 100 by 1.
j = 1;                              %j = one
for ii = 1:length(GPR_mc.time)-1        %from one through length of the GPR_mc.time - 1.  %%Why this iteration??
    if GPR_mc.time(ii+1)-GPR_mc.time(ii)>1  %After each trial there is a pause of 1 second that is used to find the indexes
        A = GPR_mc.time(ii);            %
        i_GPR_mc(j) = find(GPR_mc.time==A); %
        j=j+1;                      %
    end    
end
i_GPR_mc = nonzeros(i_GPR_mc);              %get the nonzero elements in i_GPR_mc, put them as i_GPR_mc

% 47 indexes are expected and are checked. VI SKAL TJEKKE HVOR 
if length(i_GPR_mc) == 31               %if the expected length of the array i_GPR_mc is 47          
    disp('Yes')                     %display yes
else                                %if the expected length of the array i_GPR_mc, is NOT 47, 
    disp('No, there is a mistake')  %display mistake message
end
i_GPR_mc = i_GPR_mc+1; i_GPR_mc = [i_GPR_mc;length(GPR_mc.time)]; %Add final index

% Show performance: Time it took to reach each target
P = zeros(1,32);                    %create an array of all zeros with dim 1 by 48 called P
for jj = 1:length(P)                %for loop with jj going from 1 through the length of P  %%why is this?? 
    if jj == 1                      %when jj equals 1               %%explain this loop please
        P(jj) = GPR_mc.time(i_GPR_mc(jj)-1)-1;  
    elseif jj == 32                 
        P(jj) = GPR_mc.time(i_GPR_mc(32))-GPR_mc.time(i_GPR_mc(jj-1))-1;  
    else
        P(jj) = GPR_mc.time(i_GPR_mc(jj)-1)-GPR_mc.time(i_GPR_mc(jj-1))-1;
    end
end
Score_GPR_mc = sum(P);          %Summation of all times is the performance score. 
                            %Smaller number = better (faster time at reaching the target)

% GPR (Same is done as in GPR_mc)
i_gpr_uc = zeros(100,1);
j = 1;
for ii = 1:length(GPR_uc.time)-1
    if GPR_uc.time(ii+1)-GPR_uc.time(ii)>1
        A = GPR_uc.time(ii);
        i_gpr_uc(j) = find(GPR_uc.time==A);
        j=j+1;
    end
end
i_gpr_uc = nonzeros(i_gpr_uc);

if length(i_gpr_uc) == 31
    disp('Yes')
else 
    disp('No, there is a mistake')
end
i_gpr_uc = i_gpr_uc+1; i_gpr_uc  = [i_gpr_uc;length(GPR_uc.time)];

% Show performance: Time it took to reach each target
P = zeros(1,32);
for jj = 1:length(P)
    if jj == 1
        P(jj) = GPR_uc.time(i_gpr_uc(jj)-1)-1;
    elseif jj == 32
        P(jj) = GPR_uc.time(i_gpr_uc(32))-GPR_uc.time(i_gpr_uc(jj-1))-1;
    else
        P(jj) = GPR_uc.time(i_gpr_uc(jj)-1)-GPR_uc.time(i_gpr_uc(jj-1))-1;
    end
end
Score_GPR_uc = sum(P);     %Summation of all times is the performance score. 
                        %Smaller number = better (faster time at reaching the target)

% Percentage of GPR being relatively faster. When negative, GPR_mc was faster.
B = (Score_GPR_mc - Score_GPR_uc)/Score_GPR_mc;       

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
GPR_uc(1).path.opt = zeros(length(i_gpr_uc),1);       %create an array of zeros with dim length(i_gpr) by 1 called GPR(1).path.opt
GPR_mc(1).path.opt = zeros(length(i_GPR_mc),1);         %create an array of zeros with dim length(i_GPR_mc) by 1 called GPR_mc(1).path.opt
GPR_uc(1).path.taken = zeros(length(i_gpr_uc),1);     %create an array of zeros with dim length(i_gpr) by 1 called GPR(1).path.taken
GPR_mc(1).path.taken = zeros(length(i_GPR_mc),1);       %create an array of zeros with dim length(i_GPR_mcr) by 1 called GPR_mcR(1).path.taken
for jj = 1:32   %run from jj equal to 1 through 48.    Look into the equations below for Path optimising
    GPR_uc.path.opt(jj) = sqrt(GPR_uc.target(i_gpr_uc(jj)-1,1)^2+GPR_uc.target(i_gpr_uc(jj)-1,2)^2)-0.075; %0.075 is the radius of the target
    GPR_mc.path.opt(jj) = sqrt(GPR_mc.target(i_GPR_mc(jj)-1,1)^2+GPR_mc.target(i_GPR_mc(jj)-1,2)^2)-0.075;
    if jj == 1  %if jj equals to 1, do:
        gprOut_uc = [0,0;GPR_uc.sysOut(1:i_gpr_uc(jj)-1,1),GPR_uc.sysOut(1:i_gpr_uc(jj)-1,2)];     %GPR out is set to be concatenate of 0,0 and GPR.sysOut(1:i_gpr(jj)-1,1),GPR.sysOut(1:i_gpr(jj)-1,2)
        GPROut_mc = [0,0;GPR_mc.sysOut(1:i_GPR_mc(jj)-1,1),GPR_mc.sysOut(1:i_GPR_mc(jj)-1,2)];          %GPR_mc out is set to be concatenate of 0,0 and GPR_mc.sysOut(1:i_GPR_mc(jj)-1,1),GPR_mc.sysOut(1:i_GPR_mc(jj)-1,2)
        tmax_gpr_uc = length(gprOut_uc);          %returns the length of the longest array dimension in gprOut
        tmax_GPR_mc = length(GPROut_mc);            %returns the length of the longest array dimension in GPROut_mc
        t_gpr_uc = [1:dt:tmax_gpr_uc];            %t_gpr = [start:step:stop] -> start at 1, go 0.01 each step, end at tmax_gpr
        t_GPR_mc = [1:dt:tmax_GPR_mc];              %t_GPR_mc = [start:step:stop] -> start at 1, go 0.01 each step, end at tmax_GPR_mc
        gprOut_uc_new = interp1(gprOut_uc, t_gpr_uc, 'pchip');   %set gprOut_new to be 1D interpolation of gprOut and t_gpr with 
                                                        %method pchip (Piecewise Cubic Hermite Interpolating Polynomial)
                                                        %connecting the dots
        GPROut_mc_new = interp1(GPROut_mc, t_GPR_mc, 'pchip');      %set GPROut_mc_new to be 1D interpolation of GPROut_mc and t_GPR_mc with 
                                                        %method pchip (Piecewise Cubic Hermite Interpolating Polynomial)
        count = 2;
        for kk = 3:length(gprOut_uc_new)           %for loop when kk equals 3 through the length of gprOut_new
            if sqrt((GPR_uc.target(i_gpr_uc(jj)-1,1)-gprOut_uc_new(end-(kk-1),1))^2+(GPR_uc.target(i_gpr_uc(jj)-1,2)-gprOut_uc_new(end-(kk-1),2))^2) < 0.075 
                %above: if the ?? is below the radius of the target. ?? dont understand
                count = count + 1;              %add 1 to count, put this as the new count
            else
                break
            end
        end
        real_end_gpr_uc = length(gprOut_uc_new)-count+1;  %real_end_gpr equals the length of gprOut_new - (count+1)
        d_gpr_uc = [gprOut_uc_new(1,1),gprOut_uc_new(1,2);diff(gprOut_uc_new(1:real_end_gpr_uc,:))]; %?? dont understand
        count = 2;                              %set count to be equal to 2
        for kk = 3:length(GPROut_mc_new)            %for loop when kk equals 3 through length of GPROut_mc_new
            if sqrt((GPR_mc.target(i_GPR_mc(jj)-1,1)-GPROut_mc_new(end-(kk-1),1))^2+(GPR_mc.target(i_GPR_mc(jj)-1,2)-GPROut_mc_new(end-(kk-1),2))^2) < 0.075
                %above: if the ?? is below the radius of the target. ?? dont understand
                count = count + 1;              %add 1 to count, put this as the new count
            else
                break
            end
        end
        real_end_GPR_mc = length(GPROut_mc_new)-count+1;    %real_end_GPR_mc equals the length of GPROut_mc_new - (count+1)
        d_GPR_mc = [GPROut_mc_new(1,1),GPROut_mc_new(1,2);diff(GPROut_mc_new(1:real_end_GPR_mc,:))]; %?? dont understand
    elseif jj == 32                 %else if jj is equal to 48 do:
        gprOut_uc = [GPR_uc.sysOut(i_gpr_uc(jj-1):i_gpr_uc(32),1),GPR_uc.sysOut(i_gpr_uc(jj-1):i_gpr_uc(32),2)];
        GPROut_mc = [GPR_mc.sysOut(i_GPR_mc(jj-1):i_GPR_mc(32),1),GPR_mc.sysOut(i_GPR_mc(jj-1):i_GPR_mc(32),2)];
        tmax_gpr_uc = length(gprOut_uc);  %set tmax_gpr to be the length of gprOut
        tmax_GPR_mc = length(GPROut_mc);    %set tmax_GPR_mc to be the length of GPROut_mc
        t_gpr_uc = [1:dt:tmax_gpr_uc];            %t_gpr = [start:step:stop] -> start at 1, go 0.01 each step, end at tmax_gpr
        t_GPR_mc = [1:dt:tmax_GPR_mc];              %t_GPR_mc = [start:step:stop] -> start at 1, go 0.01 each step, end at tmax_GPR_mc
        gprOut_uc_new = interp1(gprOut_uc, t_gpr_uc, 'pchip');   %set gprOut_new to be 1D interpolation of gprOut and t_gpr with 
                                                        %method pchip (Piecewise Cubic Hermite Interpolating Polynomial)
                                                        %connecting the dots
        GPROut_mc_new = interp1(GPROut_mc, t_GPR_mc, 'pchip');      %set GPROut_mc_new to be 1D interpolation of GPROut_mc and t_GPR_mc with 
                                                        %method pchip (Piecewise Cubic Hermite Interpolating Polynomial)
        count = 2;
        for kk = 3:length(gprOut_uc_new)
            if sqrt((GPR_uc.target(i_gpr_uc(jj)-1,1)-gprOut_uc_new(end-(kk-1),1))^2+(GPR_uc.target(i_gpr_uc(jj)-1,2)-gprOut_uc_new(end-(kk-1),2))^2) < 0.075
                count = count + 1;
            else
                break
            end
        end
        real_end_gpr_uc = length(gprOut_uc_new)-count+1;
        d_gpr_uc = [gprOut_uc_new(1,1),gprOut_uc_new(1,2);diff(gprOut_uc_new(1:real_end_gpr_uc,:))];
        count = 2;
        for kk = 3:length(GPROut_mc_new)
            if sqrt((GPR_mc.target(i_GPR_mc(jj)-1,1)-GPROut_mc_new(end-(kk-1),1))^2+(GPR_mc.target(i_GPR_mc(jj)-1,2)-GPROut_mc_new(end-(kk-1),2))^2) < 0.075
                count = count + 1;
            else
                break
            end
        end
        real_end_GPR_mc = length(GPROut_mc_new)-count+1;
        d_GPR_mc = [GPROut_mc_new(1,1),GPROut_mc_new(1,2);diff(GPROut_mc_new(1:real_end_GPR_mc,:))];
    else
        gprOut_uc = [GPR_uc.sysOut(i_gpr_uc(jj-1):i_gpr_uc(jj)-1,1),GPR_uc.sysOut(i_gpr_uc(jj-1):i_gpr_uc(jj)-1,2)];
        GPROut_mc = [GPR_mc.sysOut(i_GPR_mc(jj-1):i_GPR_mc(jj)-1,1),GPR_mc.sysOut(i_GPR_mc(jj-1):i_GPR_mc(jj)-1,2)];
        tmax_gpr_uc = length(gprOut_uc);
        tmax_GPR_mc = length(GPROut_mc);
        t_gpr_uc = [1:dt:tmax_gpr_uc];
        t_GPR_mc = [1:dt:tmax_GPR_mc];
        gprOut_uc_new = interp1(gprOut_uc, t_gpr_uc, 'pchip');
        GPROut_mc_new = interp1(GPROut_mc, t_GPR_mc, 'pchip');
        count = 2;
        for kk = 3:length(gprOut_uc_new)
            if sqrt((GPR_uc.target(i_gpr_uc(jj)-1,1)-gprOut_uc_new(end-(kk-1),1))^2+(GPR_uc.target(i_gpr_uc(jj)-1,2)-gprOut_uc_new(end-(kk-1),2))^2) < 0.075
                count = count + 1;
            else
                break
            end
        end
        real_end_gpr_uc = length(gprOut_uc_new)-count+1;
        d_gpr_uc = [gprOut_uc_new(1,1),gprOut_uc_new(1,2);diff(gprOut_uc_new(1:real_end_gpr_uc,:))];
        count = 2;
        for kk = 3:length(GPROut_mc_new)
            if sqrt((GPR_mc.target(i_GPR_mc(jj)-1,1)-GPROut_mc_new(end-(kk-1),1))^2+(GPR_mc.target(i_GPR_mc(jj)-1,2)-GPROut_mc_new(end-(kk-1),2))^2) < 0.075
                count = count + 1;
            else
                break
            end
        end
        real_end_GPR_mc = length(GPROut_mc_new);
        d_GPR_mc = [GPROut_mc_new(1,1),GPROut_mc_new(1,2);diff(GPROut_mc_new(1:real_end_GPR_mc,:))];
    end
    GPR_uc.path.taken(jj) = sum(sqrt(sum(d_gpr_uc.*d_gpr_uc,2)));
    GPR_mc.path.taken(jj) = sum(sqrt(sum(d_GPR_mc.*d_GPR_mc,2)));
end
clear count d_gpr d_GPR_mc dt gprOut gprOut_new jj kk GPROut_mc GPROut_mc_new real_end_gpr real_end_GPR_mc t_gpr t_GPR_mc tmax_gpr tmax_GPR_mc
% Path efficiency calculated
path_ratio_gpr_uc=(GPR_uc.path.taken - GPR_uc.path.opt)./GPR_uc.path.opt;
path_ratio_GPR_mc=(GPR_mc.path.taken - GPR_mc.path.opt)./GPR_mc.path.opt;

%% Diagonal and non-diagonal path efficiencies
% Find the indexes of diagonal and non-diagonal
% GPR
gpr_targets_uc = GPR_uc.target([1;i_gpr_uc(1:end-1)],:);
gpr_targets_check_uc = zeros(length(gpr_targets_uc),1);
for jj = 1:length(gpr_targets_uc)
    if abs(gpr_targets_uc(jj,1)) == abs(gpr_targets_uc(jj,2))
        gpr_targets_check_uc(jj) = 1;
    end
end
i_gpr_d_uc = find(gpr_targets_check_uc);          %find the indexes and values of the nonzero elements in gpr_targets_check
i_gpr_nd_uc = find(gpr_targets_check_uc==0);      %find the index equal to 0 in the gpr_targets_check

% GPR_mc
GPR_mc_targets = GPR_mc.target([1;i_GPR_mc(1:end-1)],:);
GPR_mc_targets_check = zeros(length(GPR_mc_targets),1);
for jj = 1:length(GPR_mc_targets)
    if abs(GPR_mc_targets(jj,1)) == abs(GPR_mc_targets(jj,2))
        GPR_mc_targets_check(jj) = 1;
    end
end
i_GPR_mc_d = find(GPR_mc_targets_check);            %find the indexes and values of the nonzero elements in GPR_mc_targets_check
i_GPR_mc_nd = find(GPR_mc_targets_check==0);        %find the index equal to 0 in the GPR_mc_targets_check

% Path ratios of diagonal and non-diagonal targets for GPR_mc and gpr.
path_ratio_gpr_d_uc=(GPR_uc.path.taken(i_gpr_d_uc) - GPR_uc.path.opt(i_gpr_d_uc))./GPR_uc.path.opt(i_gpr_d_uc);
path_ratio_gpr_nd_uc=(GPR_uc.path.taken(i_gpr_nd_uc) - GPR_uc.path.opt(i_gpr_nd_uc))./GPR_uc.path.opt(i_gpr_nd_uc);
path_ratio_GPR_mc_d=(GPR_mc.path.taken(i_GPR_mc_d) - GPR_mc.path.opt(i_GPR_mc_d))./GPR_mc.path.opt(i_GPR_mc_d);
path_ratio_GPR_mc_nd=(GPR_mc.path.taken(i_GPR_mc_nd) - GPR_mc.path.opt(i_GPR_mc_nd))./GPR_mc.path.opt(i_GPR_mc_nd);

%% Performance matrices: total and diagonal vs non-diagonal
% All performance measures are put into 1 matrix. This is done manually for
% each subject.

%Performance = zeros(10,6);
Performance(subject,1) = Score_GPR_mc;
Performance(subject,2) = Score_GPR_uc;
Performance(subject,3) = B;
Performance(subject,4) = mean(path_ratio_GPR_mc);
Performance(subject,5) = mean(path_ratio_gpr_uc);
Performance(subject,6) = (Performance(subject,4)-Performance(subject,5))./Performance(subject,4);

%Performance_d = zeros(10,6);
Performance_d(subject,1) = mean(path_ratio_gpr_d_uc);
Performance_d(subject,2) = mean(path_ratio_gpr_nd_uc);
Performance_d(subject,3) = mean(path_ratio_GPR_mc_d);
Performance_d(subject,4) = mean(path_ratio_GPR_mc_nd);
Performance_d(subject,5) = (Performance_d(subject,3)-Performance(subject,1))./Performance(subject,3);
Performance_d(subject,6) = (Performance_d(subject,4)-Performance(subject,2))./Performance(subject,4);

%% EXTRAS: This include trajectory simultion and a sample time check
% % Simulate Trajectory 
% close all
% hfig = figure;
% set(hfig, 'units', 'normalized','position', [0 0 1 1],'menubar', 'none','renderer','painters','name','Experiment','numbertitle','off')
% % Add target marker
% htarget = plot(GPR_mc.target(i_GPR_mc(jj)-1,1), GPR_mc.target(i_GPR_mc(jj)-1,2), 'ro','markersize',50,'linewidth', 2, 'buttondownfcn', 'uiresume'); hold on
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
% x1 = GPR_mc.sysOut(i_GPR_mc(jj-1):i_GPR_mc(jj)-1,1);
% y1 = GPR_mc.sysOut(i_GPR_mc(jj-1):i_GPR_mc(jj)-1,2);
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
% dt_GPR_mc = diff(GPR_mc.time(i(1):i(2)));diff(GPR_mc.time(i(2)+1:i(3)));
%diff(GPR_mc.time(i(3)+1:i(4)));diff(GPR_mc.time(i(4)+1:i(5)));
%diff(GPR_mc.time(i(5)+1:i(6)));diff(GPR_mc.time(i(6)+1:i(7)));
%diff(GPR_mc.time(i(7)+1:i(8)));diff(GPR_mc.time(i(8)+1:i(9)));
%diff(GPR_mc.time(i(9)+1:i(10)));diff(GPR_mc.time(i(10)+1:i(11)));
%diff(GPR_mc.time(i(11)+1:i(12)));diff(GPR_mc.time(i(12)+1:i(13)));
%diff(GPR_mc.time(i(13)+1:i(14)));diff(GPR_mc.time(i(14)+1:i(15)));
%diff(GPR_mc.time(i(15)+1:i(16)));diff(GPR_mc.time(i(16)+1:i(17)));
%diff(GPR_mc.time(i(17)+1:i(18)));diff(GPR_mc.time(i(18)+1:i(19)));diff(GPR_mc.time(i(19)+1:i(20)));
% 
% disp('The mean dt for GPR and GPR_mc are: ',mean(dt_gpr),' ',mean(dt_GPR_mc));
% disp('The std of the dt for GPR and GPR_mc are: ',std(dt_gpr),' ',std(dt_GPR_mc));

