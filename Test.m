% clc; clear; close all;
%% The Fitts law test - by Tom Baumeister
%
% This is the script used for the real experiment. The script was runned
% twice for each subject: Ones for GPR and ones for LR. Make sure to have a
% folder for each subject. Data will be saved into that folder. The data
% consists of the cursor location, target location, predicted target
% values, time and p_values.

%% Set Myoband ready to collect EMG data
  % The following 2 commented lines of codes builds the myomex file.
  % Uncommend and run if necessary.

%  sdk_path = 'C:\myo-sdk-win-0.9.0'; % root path to Myo SDK
%  build_myo_mex(sdk_path); % builds myo_mex
 
if exist('mm')
    mm.delete
end
countMyos = 1; % 1 myoband
mm = MyoMex(countMyos);
pause(1);
m1 = mm.myoData(1); % Object that collects data from the Myoband
m1.stopStreaming(); % Streaming is stopped
m1.clearLogs(); % The log file that saves the collected data is cleared.

% Load regression models
disp('Open the regression models to be used.')
if ~exist('gprMdl_dof1') || ~exist('LRmdl_1') %Skal indkommenteres når vi har linær model. 
    uiopen('*.mat')
end

% Load txt file that determines target positions and size
fileID = fopen('Fitt_real_ny.txt','r');            %load the file from Fitt_real.txt
Target = fscanf(fileID,'%f %f %f %f',[4 inf])'; %fscanf reads data from the fileID, 
                                                %with format %f = floating-point numbers. 
                                                %Reads 4 rows and to the end of the file in rows
Target = [Target;Target];                       % Double the length of the task
Target = Target(randperm(length(Target)),:);    % Randomize it


%% Initialize important parameters
Reach_time = 15; % Max time to reach target [s]
Time = Reach_time*length(Target); % Max duration of simulation [s]
samples_win = 20; % Number of samples used for prediction after windowing
dt = samples_win/200; % Sample time after windowing [s]
threshold = 0.1; % Threshold
W = Target(1,3); % Width or difficulty of the target
Dwell = 1; % Dwell time [s]
N = Time/dt; % Number of iterations

% Choose algorithm
algo = input('What algorithm? LR/GPR [GPR]: ','s');
if ~strcmp(algo,"LR")
    algo = "GPR";
    dataB(1).p = zeros(Time/dt,4); % For GPR p-values are saved
end

%% Initialization
% Matrices/vectors
window1_data = zeros(samples_win,8);    %window1 data array with only zeros of size (samples_win,8)
window2_data = zeros(samples_win,8);    %window2 data array with only zeros of size (samples_win,8)
window3_data = zeros(samples_win,8);    %window3 data array with only zeros of size (samples_win,8)
rms_data = zeros(Time/dt,8);            %rootmeansquared data array with only zeros of size (time / dt, 8)
dataB(1).time = zeros(Time/dt,1);        %the time value from data(1) is put in array with only zeros of size (time / dt, 1)
dataB(1).dof = zeros(Time/dt,4);         %the dof value from data(1) is put in array with only zeros of size (time / dt, 4)
dataB(1).target = zeros(Time/dt,2);      %the target value from data(1) is put in array with only zeros of size (time / dt, 2)
dataB(1).cursor = zeros(Time/dt,2);      %the cursor value from data(1) is put in array with only zeros of size (time / dt, 2)
dataB(1).sysOut = zeros(Time/dt,2);      %the sysOut value from data(1) is put in array with only zeros of size (time / dt, 2)

% Counters
counter_target = 1;                 %set the counter for target to be 1
counter_dwell = 0;                  %set the counter for dwell to be 0
counter_score = 0;                  %set the counter for score to be 0
counter_fail = 0;                   %set the counter for fail to be 0
counter_reach = 0;                  %set the counter for reach to be 0

% Control system
control_type = input('What control? pos/vel [pos]: ','s');
if ~strcmp(control_type,'vel')      %string compare controltype velocity
    control_type = "pos";           %with controltype position
end                                 %end comparison
if control_type == "pos"            %if controltype is position
    sys = tf(1, 1);                 %set system to be the transferfunction for position control
else                                %if controltype is not pos
    sys = tf(0.8, [0 1 0]);         %set system to be the transferfunction for velocity control
end
% State space model of transfer function
sysDisc = c2d(ss(sys),dt);          %convert statespacd model from continuous time to discrete
sysA = sysDisc.A;                   %sysA is set to be the system discrete at A
sysB = sysDisc.B;                   %sysB is set to be the system discrete at B
sysC = sysDisc.C;                   %sysC is set to be the system discrete at C
sysD = sysDisc.D;                   %sysD is set to be the system discrete at D
sysState = zeros(size(sysB));       %Creates an array of only zeros with the size of sysB, called sysState.

%% Calculate the confidence baseline 
% tilføjet af Steff
baseCI_1 = mean(o1(:,2))-mean(o1(:,1)); %Den gennemsnitlige forskel mellem nedre og øvre CI
baseCI_2 = mean(o2(:,2))-mean(o2(:,1));
baseCI_3 = mean(o3(:,2))-mean(o3(:,1));
baseCI_4 = mean(o4(:,2))-mean(o4(:,1));
baseCI_5 = mean(o5(:,2))-mean(o5(:,1));
baseCI_6 = mean(o6(:,2))-mean(o6(:,1));

%% Initialize fitts law task
close all
hfig = figure; %make a figure called hfig       [ 0 0 1 1] = [x y width height] = full screen
set(hfig, 'units', 'normalized','position', [0 0 1 1],'menubar', 'none','renderer','painters','name','Experiment','numbertitle','off')

%Barplots of the confidence - tilføjet af Mat
subplotConfi = subplot(1,2,1, 'position', [0.02 0.02 0.22 0.95]);
axis([-1 1 -1 1])
set(subplotConfi,'color', [.98,.98,.98], 'xtick', 0, 'ytick', 0);

htext_extFlex = text(-0.78, 0.90, {'Extension          Flexion'},'fontsize', 16, 'color', 'k'); 
con_ext = rectangle('Position',[-0.70 0.35 0.5 0.5],'EdgeColor',[.4,.4,.4],'Linewidth',1.2); hold on;
con_flex = rectangle('Position',[0.20 0.35 0.5 0.5],'EdgeColor',[.4,.4,.4],'Linewidth',1.2); hold on;

htext_rdUd = text(-0.78, 0.25, {'Radial Dev.    Ulnar Dev.'},'fontsize', 16, 'color', 'k'); 
con_rd = rectangle('Position',[-0.70 -0.30 0.5 0.5],'EdgeColor',[.4,.4,.4],'Linewidth',1.2); hold on;
con_ud = rectangle('Position',[0.20 -0.30 0.5 0.5],'EdgeColor',[.4,.4,.4],'Linewidth',1.2); hold on;

htext_supPro = text(-0.78, -0.40, {'Supination       Pronation'},'fontsize', 16, 'color', 'k'); 
con_sup= rectangle('Position',[-0.70 -0.95 0.5 0.5],'EdgeColor',[.4,.4,.4],'Linewidth',1.2);
con_pro = rectangle('Position',[0.20 -0.95 0.5 0.5],'EdgeColor',[.4,.4,.4],'Linewidth',1.2); hold on;

%Plotting the confidence (fill) - tilføjet af Mat
%The fill will represent the confidence by varying the height of a green

ext_fill = 0; flex_fill= 0; rd_fill= 0; ud_fill= 0; sup_fill=0; pro_fill=0;
con_ext_fill = rectangle('Position',[-0.70 0.35 0.5 ext_fill],'FaceColor','green','EdgeColor',[.4,.4,.4],'Linewidth',1.2); hold on;
con_flex_fill = rectangle('Position',[0.20 0.35 0.5 flex_fill],'FaceColor','green','EdgeColor',[.4,.4,.4],'Linewidth',1.2); hold on;

con_rd_fill = rectangle('Position',[-0.70 -0.30 0.5 rd_fill],'FaceColor','green','EdgeColor',[.4,.4,.4],'Linewidth',1.2); hold on;
con_ud_fill = rectangle('Position',[0.20 -0.30 0.5 ud_fill],'FaceColor','green','EdgeColor',[.4,.4,.4],'Linewidth',1.2); hold on;

con_sup_fill = rectangle('Position',[-0.70 -0.95 0.5 sup_fill],'FaceColor','green','EdgeColor',[.4,.4,.4],'Linewidth',1.2);
con_pro_fill = rectangle('Position',[0.20 -0.95 0.5 pro_fill],'FaceColor','green','EdgeColor',[.4,.4,.4],'Linewidth',1.2); hold on;

%Supplot for fitts law
subplotFitts = subplot(1,2,2, 'position', [0.25 0.02 0.65 0.95]);
% Add target marker and cross
htarget = plot(Target(1,1), Target(1,2), 'ro','markersize',5,'MarkerFaceColor', [0.9100 0.4100 0.1700], 'linewidth', 2, 'buttondownfcn', 'uiresume'); hold on
htarget_cross = plot(Target(1,1), Target(1,2), 'k+','markersize',Target(1,4),'linewidth', .7, 'buttondownfcn', 'uiresume'); hold on
% Add position cursor (of controlled system)
hpos = plot(0, 0, 'ko', 'markersize', 12, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0.1 0.1 0.1], 'linewidth', 2, 'buttondownfcn', 'uiresume');
axis([-1 1 -1 1])
set(subplotFitts, 'color', [.98,.98,.98],'DataAspectRatio',[1 1 1], 'units', 'normalized', 'position', [0.25 0.02 0.65 0.95], 'xtick', 0, 'ytick', 0, 'LineWidth', 3,'GridColor','k')%,'drawmode', 'fast')
grid on
% Some information shown during the task
htext = text(0, 0.7, {'Move the black circular cursor to the orange target';''; 'Click the cursor to start...'}, 'fontsize', 22, 'color', 'k', 'horizontalalignment', 'center');
htext_score1 = text(1.01, 0.8, {'Score:'},'fontsize', 14, 'color', 'k'); htext_score2 = text(1.2, 0.8, num2str(counter_score),'fontsize', 14, 'color', 'k');
htext_score3 = text(1.01, 0.7, {'Time:'},'fontsize', 14, 'color', 'k'); htext_score4 = text(1.2, 0.7, num2str('-'),'fontsize', 14, 'color', 'k');
htext_score5 = text(1.01, 0.6, {'Fails:'},'fontsize', 14, 'color', 'k'); htext_score6 = text(1.2, 0.6, num2str(counter_fail),'fontsize', 14, 'color', 'k');



% Wait for start
uiwait
% Hide mouse pointer
set(gcf, 'pointer', 'custom','pointershapecdata', nan(16,16))
delay = 3; % Countdown delay [s]
% Countdown
set(htext, 'string', num2str(delay));
tic        %start stopwatch clock
for ll = 1:delay
    drawnow
    while (toc < ll)
    end
    set(htext,'string', num2str(delay-ll))
end
set(htext, 'string', '')
drawnow



%% START SIMULATION
kk = samples_win;               %set the samples in window to a variable called kk
t_score = tic;                  %set the stopwatch values to a variable called t_score (to show the subject the time during test)
tic
m1.startStreaming();            %switch the streammode to start for the myo armband 1 (m1)
for ii = 1:N
    % Wait until log is filled with enough samples
    while length(m1.emg_log)<kk %while the length of emg log for m1 is lower than the kk
    end
    m1.stopStreaming();         %switch the streammode to stop for m1
    
    % Derive Rootmeansquared (RMS) value of window datas
    window3_data(1:samples_win,:) = m1.emg_log(kk-(samples_win-1):kk,:);
    rms_data(N,:)=rms([ window1_data(1:samples_win,:) ; window2_data(1:samples_win,:) ; window3_data(1:samples_win,:) ]);
    
    % Target predictions: GPR or LR
    switch algo
        case "GPR"              %Gaussian process regression
            [dof1,p1,ci1] = predict(gprMdl_dof1,rms_data(N,:));
            [dof2,p2,ci2] = predict(gprMdl_dof2,rms_data(N,:));
            [dof3,p3,ci3] = predict(gprMdl_dof3,rms_data(N,:));
            [dof4,p4,ci4] = predict(gprMdl_dof4,rms_data(N,:));
            [dof5,p5,ci5] = predict(gprMdl_dof5,rms_data(N,:));
            [dof6,p6,ci6] = predict(gprMdl_dof6,rms_data(N,:));
            dataB.p(ii,1:6)=[p1,p2,p3,p4,p5,p6];
            
            
        case "LR"               %Linear regression 
              dof1 = predict(LRmdl_1,rms_data(N,:));
              dof2 = predict(LRmdl_2,rms_data(N,:));
              dof3 = predict(LRmdl_3,rms_data(N,:));
              dof4 = predict(LRmdl_4,rms_data(N,:));
              dof5 = predict(LRmdl_5,rms_data(N,:));
              dof6 = predict(LRmdl_6,rms_data(N,:));
    end
    
     % Predictions are ranging [0 1]. Set to [-1 1] for x and y
    if dof1 >= dof2                 %if dof1 is bigger than or equal to dof2
        dofA = -dof1;               %put dofA to be -dof1
    else
        dofA = dof2;                %else put dofA to be dof2
    end
    if dof4 >= dof6                 %if dof4 is bigger than or equal to dof6
        dofB = dof4;                %put dofB to be dof4                            Why is this reverse from MYO4?
    else
        dofB = -dof6;               %else dofB to be -dof4                          Why is this reverse from MYO4?
    end
    
    % Use movement thresholds for minima and maxima
    if abs(dofB) < threshold            %if the absolute value of dofB is lower than the threshold
        dofB=0;                         %set dofB to be equal to 0
    end
    dofB(dofB>1)=1; dofB(dofB<-1)=-1;   %if dofB is bigger than 1 set to 1. if it is less than -1, set to -1.

    if abs(dofA) < threshold            %if the absolute value of dofA is lower than the threshold
        dofA=0;                         %set dofA to be equal to 0
    end
    dofA(dofA>1)=1; dofA(dofA<-1)=-1;   %if dofA is bigger than 1 set to 1. if it is less than -1 set to -1.
    cursor = [dofA,dofB];               %set the cursor to be the concatenation of dofB and dofA (put two matrices together to create a larger one).
    
    % Update controlled system states
    sysState = sysA * sysState + sysB * cursor;     %set the system state to be the sysA times old system state plus sysB times the cursor matrix
    sysOut = sysC * sysState + sysD * cursor;       %set the system state out to be the sysC times old system state plus sysD times the cursor matrix
   
    % Limit the play field to the border
    sysOut(sysOut>1)=1;                             %If systemstate is bigger than 1, set to 1.
    sysOut(sysOut<-1)=-1;                           %If systemstate is less than -1, set to -1.
    
    % Store data                                    %Why is (ii,:)? 
    dataB.time(ii,:) = toc;                          %store the elapsed time from the stopwatch (toc=time since tic start) in data.time array
    dataB.sysOut(ii,:) = sysOut;                     %store the sysOut controlled system state out to be stored in data.sysOut array
    dataB.target(ii,:) = Target(counter_target,1:2); %store the Target data 
    dataB.cursor(ii,:) = cursor;                     %store the cursor data
    dataB.dof(ii,:) = [dof1,dof2,dof4,dof6];         %store the degrees of freedom data
    
    % Update screen
    set(hpos, 'xdata', sysOut(1), 'ydata', sysOut(2));   %set the hpos to have x-data = sysOut(1) and y-data = sysOut(2)
    drawnow  
    
    % Update confidence bars HER ER DU I GANG MATHILDE - FIND UD AF AT
    % OPDATERE FIGUREN.
    ext_fill = min(max(0.5*(1-(((ci1(:,2)-ci1(:,1))-baseCI_1)/baseCI_1)),0),0.5);
    %ext_fill = max(0.5- ci1(:,2)-dof1,0);
    
    set(con_ext_fill, 'Position', [-0.70 0.35 0.5 ext_fill]);
    
    %flex_fill = max(0.5- ci2(:,2)-dof2,0);
    flex_fill = min(max(0.5*(1-(((ci2(:,2)-ci2(:,1))-baseCI_2)/baseCI_2)),0),0.5);
    set(con_flex_fill, 'Position', [0.20 0.35 0.5 flex_fill]);
    
    %rd_fill = max(0.5- ci4(:,2)-dof4,0);
    rd_fill = min(max(0.5*(1-(((ci4(:,2)-ci4(:,1))-baseCI_4)/baseCI_4)),0),0.5);
    set(con_rd_fill, 'Position', [-0.70 -0.30 0.5 rd_fill]);
    
    %ud_fill = max(0.5- ci6(:,2)-dof6,0);
    ud_fill = min(max(0.5*(1-(((ci6(:,2)-ci6(:,1))-baseCI_6)/baseCI_6)),0),0.5);
    set(con_ud_fill, 'Position', [0.20 -0.30 0.5 ud_fill]);
    drawnow
    
    
   
    
    
    % Update target position
    if sqrt((Target(counter_target,1)-sysOut(1))^2+(Target(counter_target,2)-sysOut(2))^2) <= W && counter_dwell < Dwell/dt  % When cursor is within width of target change the color of the target
        set(htarget,'MarkerFaceColor',[0.3 .9 0.3],'MarkerEdgeColor',[0.05 .75 0.05]);
        set(htarget_cross, 'xdata', Target(counter_target,1),'ydata', Target(counter_target,2));
        counter_dwell = counter_dwell + 1;
        counter_reach = counter_reach+1;
    elseif counter_dwell >= Dwell/dt % HIT, when dwell time is reached
        counter_score=counter_score+1;
        set(htext_score2,'string',num2str(counter_score));  %update score
        set(htext_score4,'string',num2str(toc(t_score)));   %show time
        counter_dwell = 0;
        counter_reach = 0;
        counter_target=counter_target+1;        %take next target
        if counter_target >= length(Target)+1   %if counter_target is greater than or equal to the length of Target +1
            break;                              %break
        end
        W = Target(counter_target,3);           %update position of the target
        set(htarget, 'xdata', Target(counter_target,1), 'ydata', Target(counter_target,2), 'markersize', Target(counter_target,4),'MarkerFaceColor', [0.9100 0.4100 0.1700], 'MarkerEdgeColor','r');
        set(htarget_cross, 'xdata', Target(counter_target,1),'ydata', Target(counter_target,2),'markersize', Target(counter_target,4));
        % Place cursor back to position for control (only works for velocity control)
        if control_type == 'vel'
            set(hpos, 'xdata', 0, 'ydata', 0,'MarkerFaceColor', [0.6 0.6 0.6],'MarkerEdgeColor', [0.6 0.6 0.6]);
            pause(1)
            sysState = [0 0]; %is done in order to not use predictions during the 1 second wait
            t_score = tic;
            set(hpos, 'xdata', 0, 'ydata', 0,'MarkerFaceColor', [0.1 0.1 0.1],'MarkerEdgeColor', 'k');
        end
    elseif counter_reach >= Reach_time/dt   %FAIL, when maximum time has passed. For the rest it's the same code as score.
        counter_fail=counter_fail+1;        %update fail score
        set(htext_score4,'string',num2str(Reach_time));     %show the time it takes to reach
        set(htext_score6,'string',num2str(counter_fail));   %update the fails
        counter_target=counter_target+1;    %take next target
        counter_dwell = 0;                  %set the counter for dwell to be 0
        counter_reach = 0;                  %set the counter for reach to be 0
        if counter_target >= length(Target)+1   %if the target counter is greater or equal to the length of the target +1
            break;                              %break
        end
        W = Target(counter_target,3);
        set(htarget, 'xdata', Target(counter_target,1), 'ydata', Target(counter_target,2), 'markersize', Target(counter_target,4),'MarkerFaceColor', [0.9100 0.4100 0.1700], 'MarkerEdgeColor','r');
        set(htarget_cross, 'xdata', Target(counter_target,1),'ydata', Target(counter_target,2),'markersize', Target(counter_target,4));
        % Place cursor back to position for velocity control
        if control_type == 'vel'
            set(hpos, 'xdata', 0, 'ydata', 0,'MarkerFaceColor', [0.6 0.6 0.6],'MarkerEdgeColor', [0.6 0.6 0.6]);
            pause(1)
            sysState = [0 0];
            t_score = tic;
            set(hpos, 'xdata', 0, 'ydata', 0,'MarkerFaceColor', [0.1 0.1 0.1],'MarkerEdgeColor', 'k');
        end
    else
        set(htarget, 'xdata', Target(counter_target,1), 'ydata', Target(counter_target,2),'MarkerFaceColor', [0.9100 0.4100 0.1700], 'MarkerEdgeColor','r');
        set(htarget_cross, 'xdata', Target(counter_target,1),'ydata', Target(counter_target,2));
        counter_dwell = 0;
        counter_reach = counter_reach+1;
    end
    
    % For next iteration, update windows 1 and 2. Window 1 is filled in at iteration 3.
    window2_data(1:samples_win,:) = m1.emg_log(kk-(samples_win-1):kk,:);
    if ii >= 3
        window1_data(1:samples_win,:) = m1.emg_log((kk-(2*samples_win-1)):kk-samples_win,:);
    end
    kk=kk+samples_win;
    m1.startStreaming();            %change the streammode to start
end
m1.stopStreaming();                 %change the streammode to stop, to stop streaming from m1
m1.clearLogs();                     %clear the logs for m1

% END SIMULATION
set(htext, 'string', 'Done...')     %Write the text "done" when simulation is done
pause(3);                           %Pause for 3 seconds
close(hfig)                         %close figure

%% Calculate performance and save data
dataB.time = nonzeros(dataB.time);                        %returns the nonzero elements in data.time
dataB.sysOut = dataB.sysOut(1:length(dataB.time),:);       %returns the ??
dataB.target = dataB.target(1:length(dataB.time),:);       %??
dataB.cursor = dataB.cursor(1:length(dataB.time),:);       %??
dataB.dof = dataB.dof(1:length(dataB.time),:);             %??
if algo == "GPR"                                        %if algorithm case is Gaussian Process Regression
    dataB.p = dataB.p(1:length(dataB.time),:);             %save the data ??
    save('dataGPR','-struct','data')                    %save the data of GRP in a structure field in a file called dataGPR
else 
    save('dataLR','-struct','data')                     %save the data of LR in a structure field in a file called dataLR
end

%% Plots
%Plot 1 - Trajectory
hfig_traj = figure;
set(hfig_traj, 'units', 'normalized','position', [0 0 1 1],'menubar', 'none','renderer','painters','name','Experiment','numbertitle','off','Color','w')
plot(dataB.target(1:length(dataB.target),1),dataB.target(1:length(dataB.target),2),'go')
hold on
plot(dataB.sysOut(1:length(dataB.target),1),dataB.sysOut(1:length(dataB.target),2),'bx')
axis([-1 1 -1 1])
set(gca, 'color', [.98,.98,.98],'DataAspectRatio',[1 1 1], 'units', 'normalized', 'position', [0 0 1 1], 'xtick', 0, 'ytick', 0, 'LineWidth', 3,'GridColor','k')%,'drawmode', 'fast')
grid on
if algo == "LR"
    export_fig(hfig_traj,'-pdf','filename','LR_trajectory');
else
    export_fig(hfig_traj,'-pdf','filename','GPR_trajectory');
end

%Plot 2 - Predicted target values and the uncertainty values (GPR only)
%over time.
if algo == "LR"
    subplot = @(m,n,p) subtightplot(m, n, p, [0.06 0.03], [0.08 0.05], [0.08 0.03]);
    fig_plots=figure('DefaultAxesPosition', [0.1, 0.1, 0.8, 0.8], 'units', 'normalized','position', [0 0 0.8 1],'Color','w');
    ax1 = subplot(2,1,1);
    plot(1:length(dataB.cursor),dataB.cursor,'LineWidth',1.5)
    ylabel('normalized [-]')
    ylim([-1.1 1.1])
    legend('x','y')
    title('','Fontsize',14)
    ax2 = subplot(2,1,2);
    plot(1:length(dataB.dof),dataB.dof,'LineWidth',1.5)
    ylabel('normalized [-]')
    ylim([0 1.1])
    legend('dof1','dof2','dof4','dof6')
    xlabel('samples [#]')
    set([ax1,ax2],'fontsize',11)
    export_fig(fig_plots,'-pdf','filename','LR_plots');
else
    subplot = @(m,n,p) subtightplot(m, n, p, [0.06 0.03], [0.08 0.05], [0.08 0.03]);
    fig_plots=figure('DefaultAxesPosition', [0.1, 0.1, 0.8, 0.8], 'units', 'normalized','position', [0 0 0.8 1],'Color','w');
    ax1 = subplot(3,1,1);
    plot(1:length(dataB.cursor),dataB.cursor,'LineWidth',1.5)
    ylabel('normalized [-]')
    ylim([-1.1 1.1])
    legend('x','y')
    title('','Fontsize',14)
    ax2 = subplot(3,1,2);
    plot(1:length(dataB.dof),dataB.dof,'LineWidth',1.5)
    ylabel('normalized [-]')
    ylim([0 1.1])
    legend('dof1','dof2','dof4','dof6')
    ax3 = subplot(3,1,3);
    plot(1:length(dataB.p),dataB.p,'LineWidth',1.5)
    ylabel('uncertainty [ratio]')
    legend('p1','p2','p3','p4')
    xlabel('samples [#]')
    set([ax1,ax2,ax3],'fontsize',11)
    export_fig(fig_plots,'-pdf','filename','GPR_plots');
end