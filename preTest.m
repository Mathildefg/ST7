% clc; clear; 
close all;
%% A pretest script - by Tom Baumeister
%
% The script performs the fitts law test for training prior to the real test.
% This is done in order to avoid too much learning during the real test. 
% Both GPR and LR models are randomized.

%% Set Myoband ready to collect EMG data
  % The following 2 commented lines of codes builds the myomex file.
  % Uncommend and run if necessary.
  
% sdk_path = 'C:\myo-sdk-win-0.9.0'; % root path to Myo SDK
% build_myo_mex(sdk_path); % builds myo_mex
if exist('mm')
    mm.delete
end
countMyos = 1; %One myoband is used
mm = MyoMex(countMyos); %lav en ny MyoMex
pause(1);
m1 = mm.myoData(1);     %sæt m1 til at være MyoMex myoData(1). Object that collects data from the Myoband
m1.stopStreaming();     %Although we can't stop the data from being passed to in MyoData, we can
                        %toggle streaming mode by using the methods stopStreaming() and startStreaming().
m1.clearLogs();         %clear the m1 logs.

% Load regression models if they are not in the workspace yet.
% Models should be named: gprMdl_dofx and LRmdl_x in which x is the dof
% from 1 to 4.
if ~exist('gprMdl_dof1') || ~exist('LRmdl_1')
    uiopen('*.mat')
end

% Load target file. File can be adjusted. Columns represent: [x, y, width in
% which cursor needs to be, width of the target that is shown during test].
% In this txt file target consist of 24 different positions all having the
% same width and all the positions are doubled and randomized, resulting in
% 48 random target positions.
fileID = fopen('Fitt_real_ny.txt','r');                
Target = fscanf(fileID,'%f %f %f %f',[4 inf])';     %fscanf reads data from the fileID, 
                                                    %with format %f = floating-point numbers. 
                                                    %Reads 4 rows and to the end of the file in rows
Target = [Target;Target];                           %Creates a vertical concatenating array (putting one array under the other)
Target = Target(randperm(length(Target)),:);        %randperm = random permutation of the length of Target

%% Initialize important parameters
Reach_time = 15; % Max time to reach target [s]
Time = Reach_time*length(Target); % Duration of simulation [s]
samples_win = 20; % Sample time after windowing [#samples]
dt = samples_win/200; % Sample time after windowing [s]
threshold = 0.1; % Threshold
W = Target(1,3); % Width or difficulty of the target
Dwell = 1; % Dwell time [s]
N = Time/dt; % Number of iterations
target_times = 10; % Number of targets to hit for learning
dof1 = 0; dof2 = 0; dof3 = 0; dof4 = 0; dof5 = 0; dof6 = 0; dofA = 0; dofB = 0; dofC = 0;
prevdofC = 0; nydofC = 0;
xaksen = [dof1, dof2, dof3, dof4, dof5, dof6, dofA, dofB, dofC];

% Two algorithms
algo_s = ["GPR","LR"];

%% Initialization
% Matrices/vectors
window1_data = zeros(samples_win,8); % 3 windows are used
window2_data = zeros(samples_win,8); 
window3_data = zeros(samples_win,8);
rms_data = zeros(Time/dt,8); % resulting rms data as predictors

% Counters
counter_target = 1;             %set the counter for target to be 1
counter_dwell = 0;              %set the counter for dwell to be 0
counter_score = 0;              %set the counter for score to be 0
counter_fail = 0;               %set the counter for fail to be 0
counter_reach = 0;              %set the counter for reach to be 0

% Control system: position or velocity control.
control_type = input('What control? pos/vel [pos]: ','s');
if ~strcmp(control_type,"pos") && ~strcmp(control_type,'vel')
    control_type = "pos";
end
if control_type == "pos"
    sys = tf(1, 1); %Position control
else 
    sys = tf(0.8, [0 1 0]); %Velocity control
end
sysDisc = c2d(ss(sys),dt);      %convert statepsacemodel from continuous time to discrete
sysA = sysDisc.A;               %sysA is set to be the system discrete at A
sysB = sysDisc.B;               %sysB is set to be the system discrete at B
sysC = sysDisc.C;               %sysC is set to be the system discrete at C
sysD = sysDisc.D;               %sysD is set to be the system discrete at D
sysState = zeros(size(sysB));   %Creates an array of only zeros with the size of sysB, called sysState.
%% Calculate the confidence baseline 
% tilføjet af Steff
baseCI_1 = mean(o1(:,2))-mean(o1(:,1)); %Den gennemsnitlige forskel mellem nedre og øvre CI
baseCI_2 = mean(o2(:,2))-mean(o2(:,1));
baseCI_3 = mean(o3(:,2))-mean(o3(:,1));
baseCI_4 = mean(o4(:,2))-mean(o4(:,1));
baseCI_5 = mean(o5(:,2))-mean(o5(:,1));
baseCI_6 = mean(o6(:,2))-mean(o6(:,1));


%% Build figure
close all
%plot for dof
bargraph = figure;
set(bargraph, 'units', 'normalized','position', [0 0 0.3 0.3],'menubar', 'none','renderer','painters','name','BarDoF','numbertitle','off','Color','w')
b = bar(xaksen);
ylim([-0.5 0.5]);

hfig = figure;
set(hfig, 'units', 'normalized','position', [0 0 1 1],'menubar', 'none','renderer','painters','name','Experiment','numbertitle','off')

%Barplots of the confidence - tilføjet af Mat
subplotConfi = subplot(1,2,1, 'position', [0.02 0.02 0.22 0.95]);
axis([-1 1 -1 1])
set(subplotConfi,'color', [.98,.98,.98], 'xtick', 0, 'ytick', 0);

htext_ext = text(-0.70, 0.90, {'Extension'},'fontsize', 16, 'color', 'k'); 
con_ext = rectangle('Position',[-0.70 0.35 0.5 0.5],'EdgeColor',[.4,.4,.4],'Linewidth',1.2); hold on;

htext_Flex = text(0.20, 0.90, {'Flexion'},'fontsize', 16, 'color', 'k');
con_flex = rectangle('Position',[0.20 0.35 0.5 0.5],'EdgeColor',[.4,.4,.4],'Linewidth',1.2); hold on;

htext_rd = text(-0.70, 0.25, {'Radial Dev.'},'fontsize', 16, 'color', 'k');
con_rd = rectangle('Position',[-0.70 -0.30 0.5 0.5],'EdgeColor',[.4,.4,.4],'Linewidth',1.2); hold on;

htext_Ud = text(0.20, 0.25, {'Ulnar Dev.'},'fontsize', 16, 'color', 'k');
con_ud = rectangle('Position',[0.20 -0.30 0.5 0.5],'EdgeColor',[.4,.4,.4],'Linewidth',1.2); hold on;

htext_sup = text(-0.70, -0.40, {'Open'},'fontsize', 16, 'color', 'k'); 
con_sup = rectangle('Position',[-0.70 -0.95 0.5 0.5],'EdgeColor',[.4,.4,.4],'Linewidth',1.2); hold on;

htext_pro = text(0.20, -0.40, {'Close'},'fontsize', 16, 'color', 'k'); 
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

%Subplot for fitts law
subplotFitts = subplot(1,2,2, 'position', [0.25 0.02 0.65 0.95]);

% Add target marker
htarget = plot(Target(1,1), Target(1,2), 'ro','markersize',Target(1,4),'MarkerFaceColor', [0.9100 0.4100 0.1700], 'linewidth', 2, 'buttondownfcn', 'uiresume'); hold on
% Add target cross
htarget_cross = plot(Target(1,1), Target(1,2), 'k+','markersize',Target(1,4),'linewidth', .7, 'buttondownfcn', 'uiresume'); hold on
% Add position (of controlled system) marker
hpos = plot(0, 0, 'ko', 'markersize', 12, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0.1 0.1 0.1], 'linewidth', 2, 'buttondownfcn', 'uiresume');
axis([-1 1 -1 1])
set(gca, 'color', [.98,.98,.98],'DataAspectRatio',[1 1 1], 'units', 'normalized', 'position', [0 0 1 1], 'xtick', 0, 'ytick', 0, 'LineWidth', 3,'GridColor','k')%,'drawmode', 'fast')
grid on
% Some info text
htext = text(0, 0.7, {'Move the black cursor to the orange target';''; 'Click on the cursor to start...'}, 'fontsize', 24, 'color', 'k', 'horizontalalignment', 'center');
htext_score1 = text(1.2, 0.8, {'Score:'},'fontsize', 19, 'color', 'k'); htext_score2 = text(1.4, 0.8, num2str(counter_score),'fontsize', 19, 'color', 'k');
htext_score3 = text(1.2, 0.7, {'Time:'},'fontsize', 19, 'color', 'k'); htext_score4 = text(1.4, 0.7, num2str('-'),'fontsize', 19, 'color', 'k');
htext_score5 = text(1.2, 0.6, {'Fails:'},'fontsize', 19, 'color', 'k'); htext_score6 = text(1.4, 0.6, num2str(counter_fail),'fontsize', 19, 'color', 'k');
% Wait for start. When clicked on the cursor the simulation starts
uiwait
% Hide mouse pointer
set(gcf, 'pointer', 'custom','pointershapecdata', nan(16,16))
delay = 3; % Countdown delay [s]
% Countdown
set(htext, 'string', num2str(delay));
tic     %start stopwatch timer
for ll = 1:delay                    %??
    drawnow
    while (toc < ll)
    end
    set(htext,'string', num2str(delay-ll))
end
set(htext, 'string', '')
drawnow


%% START SIMULATION
kk = samples_win;
t_score = tic; % This one is used to show the time to the subjecst during test
tic % This timer is used as data
m1.startStreaming();   %change the streammode to start for m1
for ii = 1:N
    % Wait until log is filled with kk samples. This is the factor that
    % determines how fast the simulation updates.
    while length(m1.emg_log)<kk
    end
    m1.stopStreaming();  %change the streammode to stop for m1
    
    % Derive RMS value of all 3 windows. In the first iteration windows 1
    % and 2 are set to 0 and only window 3 is filled in.
    window3_data(1:samples_win,:) = m1.emg_log(kk-(samples_win-1):kk,:);
    rms_data(N,:)=rms([ window1_data(1:samples_win,:) ; window2_data(1:samples_win,:) ; window3_data(1:samples_win,:) ]);
    
    % Target predictions: either GPR or LR
    idx = randperm(2,1);   %returns a row vector containing 1 unique integer selected randomly from 1 to 2.
    algo = algo_s(idx);    %set the algorithm to be either algo_s(1) or (2) depending on randperm.
    switch algo
        case "GPR"         %if algorithm is Gaussian Process Regression
            [dof1,p1,ci1] = predict(gprMdl_dof1,rms_data(N,:));
            [dof2,p2,ci2] = predict(gprMdl_dof2,rms_data(N,:));
            [dof3,p3,ci3] = predict(gprMdl_dof3,rms_data(N,:));
            [dof4,p4,ci4] = predict(gprMdl_dof4,rms_data(N,:));
            [dof5,p5,ci5] = predict(gprMdl_dof5,rms_data(N,:));
            [dof6,p6,ci6] = predict(gprMdl_dof6,rms_data(N,:));
            dataB.p(ii,1:6)=[p1,p2,p3,p4,p5,p6];
            
        case "LR"          %if algorithm is Linear regression
              dof1 = predict(LRmdl_1,rms_data(N,:));
              dof2 = predict(LRmdl_2,rms_data(N,:));
              dof3 = predict(LRmdl_3,rms_data(N,:));
              dof4 = predict(LRmdl_4,rms_data(N,:));
              dof5 = predict(LRmdl_5,rms_data(N,:));
              dof6 = predict(LRmdl_6,rms_data(N,:));
    end
    
    
      %% confidence interval istedet for dof1
    nyci1 = ci1(1)+ci1(2);
    nyci2 = ci2(1)+ci2(2);
    nyci3 = ci3(1)+ci3(2);
    nyci4 = ci4(1)+ci4(2);
    nyci5 = ci5(1)+ci5(2);
    nyci6 = ci6(1)+ci6(2);
     
    if  nyci1 <= nyci2 %hvis confidence interval for ext er mindre end flex
        dofA = -nyci1; %sæt dofA til at være -confidence interval for ext
    else
        dofA = nyci2; %elllers sæt til confidence interval for flex
    end
    if nyci3 <= nyci5 %hvis confidence int for rad er mindre end uln
        dofB = -nyci3; %sæt dofB til at være -con int for rad
    else
        dofB = nyci5; %ellers sæt til at conf int for uln
    end
    
    %cursor size = dofC. bliver ikke brugt lige pt til noget
    if nyci4 <= nyci6 %hvis confidence interval for 3 er mindre end 5
        dofC = -nyci4;%sæt dofC til at være -ci3
    else
        dofC = nyci6;
    end

%       % Predictions are ranging [0 1]. Set to [-1 1] for x and y
%     if dof1 >= dof2                 %if dof1 is bigger than or equal to dof2
%         dofA = -dof1;               %put dofA to be -dof1
%     else
%         dofA = dof2;                %else put dofA to be dof2
%     end
%     if dof3 >= dof4                 %if dof3 is bigger than or equal to dof4
%         dofB = dof3;                %put dofB to be dof3                            Why is this reverse from previous MYO4?
%     else
%         dofB = -dof4;               %else dofB to be -dof4                          Why is this reverse from previous MYO4?
%     end
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
        
    % Update the cursor position
    set(hpos, 'xdata', sysOut(1), 'ydata', sysOut(2));
    drawnow
    
    % Update Bargraph
    xaksen = [dof1, dof2, dof3, dof4, dof5, dof6, dofA, dofB, dofC, prevdofC, nydofC];
    b.YData = xaksen;
    
    % Update confidence bars
    
    %ext_fill = max(0.5- ci1(:,2)-dof1,0);
    ext_fill = min(max(0.5*(1-(((ci1(:,2)-ci1(:,1))-baseCI_1)/baseCI_1)),0),0.5);
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
    
    %pro_fill = max(0.5- ci3(:,2)-dof3,0);
    pro_fill = min(max(0.5*(1-(((ci3(:,2)-ci3(:,1))-baseCI_3)/baseCI_3)),0),0.5);
    set(con_pro_fill, 'Position', [0.20 -0.95 0.5 pro_fill]);
    
    %sub_fill = max(0.5-ci5(:,2)-dof5,0);
    sup_fill = min(max(0.5*(1-(((ci5(:,2)-ci5(:,1))-baseCI_5)/baseCI_5)),0),0.5);
    set(con_sup_fill, 'Position', [-0.70 -0.95 0.5 sup_fill]);
    
    
    drawnow
    
    % Update target position and keep track of the score, fails and time
    if sqrt((Target(counter_target,1)-sysOut(1))^2+(Target(counter_target,2)-sysOut(2))^2) <= W && counter_dwell < Dwell/dt % When cursor is within width of target change the color of the target
        set(htarget,'MarkerFaceColor',[0.3 .9 0.3],'MarkerEdgeColor',[0.05 .75 0.05]);
        set(htarget_cross, 'xdata', Target(counter_target,1),'ydata', Target(counter_target,2));
        counter_dwell = counter_dwell + 1;
        counter_reach = counter_reach + 1;
    elseif counter_dwell >= Dwell/dt % HIT, when dwell time is reached
        counter_score = counter_score + 1;
        set(htext_score2,'string',num2str(counter_score)); % Update score
        set(htext_score4,'string',num2str(toc(t_score))); % Show time
        counter_dwell = 0;      %set the counter for dwell to be 0
        counter_reach = 0;      %set the counter for reach to be 0
        counter_target=counter_target+1; % Take next target
        if counter_target >= length(Target)+1
            break;
        end
        if counter_target-1 >= target_times
            break
        end
        W = Target(counter_target,3); % Update position of the target
        set(htarget, 'xdata', Target(counter_target,1), 'ydata', Target(counter_target,2), 'markersize', Target(counter_target,4),'MarkerFaceColor', [0.9100 0.4100 0.1700], 'MarkerEdgeColor','r');
        set(htarget_cross, 'xdata', Target(counter_target,1),'ydata', Target(counter_target,2),'markersize', Target(counter_target,4));
        % Place cursor back to position (only works for velocity control)
        if control_type == 'vel'
            set(hpos, 'xdata', 0, 'ydata', 0,'MarkerFaceColor', [0.6 0.6 0.6],'MarkerEdgeColor', [0.6 0.6 0.6]);
            pause(1)
            sysState = [0 0]; % In order to not use predictions during the 1 second wait
            t_score = tic;
            set(hpos, 'xdata', 0, 'ydata', 0,'MarkerFaceColor', [0.1 0.1 0.1],'MarkerEdgeColor', 'k');
        end
    elseif counter_reach >= Reach_time/dt   %FAIL, when maximum time has passed. For the rest it's the same code as score.
        counter_fail = counter_fail + 1;    %Update fail score
        set(htext_score4,'string',num2str(Reach_time));     %show the time it takes to reach
        set(htext_score6,'string',num2str(counter_fail));   %update the fails
        counter_target=counter_target+1;    %take next target
        counter_dwell = 0;                  %set the counter for dwell to be 0
        counter_reach = 0;                  %set the counter for reach to be 0
        if counter_target >= length(Target)+1   %if the target counter is greater or equal to the length of the target +1
            break;                              %break
        end
        if counter_target-1 >= target_times     %if the target counter -1 is greater or equal to the target_times
            break                               %break
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
        set(htarget, 'xdata', Target(counter_target,1), 'ydata', Target(counter_target,2),'MarkerFaceColor', [0.9100 0.4100 0.1700], 'MarkerEdgeColor','r'); %When cursor is not within target, it updates color
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