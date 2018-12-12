%% Bar plot of Fitts's law test results

S = std(data);
m = mean(data);

TP = [m(1),m(2),m(3)];
PE = [m(4),m(5),m(6)];
CR = [m(7),m(8),m(9)];

TP_errorplot = [m(1);m(2);m(3)];
PE_errorplot = [m(4);m(5);m(6)];
CR_errorplot = [m(7);m(8);m(9)];

x = [1;2;3]; 
 
std_TP =[S(1);S(2);S(3)];
std_PE =[S(4);S(5);S(6)];
std_CR =[S(7);S(8);S(9)];
 
figure
hold on

subplot(1,3,1)
hold on
set(gca,'FontSize',14) % Creates an axes and sets its FontSize to 18
b = bar(TP);
b.FaceColor = 'flat';
b.CData(2,:) = [0.7 0.3 0.3];
b.CData(3,:) = [0.3 0.7 0.3];
errorbar(x,TP_errorplot,std_TP,'k.')
xticks([1 2 3])
xticklabels({'LR' 'GPR' 'GPRc'})
xlabel('Regression model','FontSize',16)
ylabel('Throughput')
%legend('Low CI and high output','High CI and low output','Low CI and low output','High CI and high output','STD')

subplot(1,3,2)
hold on
set(gca,'FontSize',14) % Creates an axes and sets its FontSize to 18
b = bar(PE);
b.FaceColor = 'flat';
b.CData(2,:) = [0.7 0.3 0.3];
b.CData(3,:) = [0.3 0.7 0.3];
errorbar(x,PE_errorplot,std_PE,'k.')
xticks([1 2 3])
xticklabels({'LR' 'GPR' 'GPRc'})
xlabel('Regression model','FontSize',16)
ylabel('Path efficiency')
title('Performance measures','FontSize',20)

subplot(1,3,3)
hold on
set(gca,'FontSize',14) % Creates an axes and sets its FontSize to 18
b = bar(CR);
b.FaceColor = 'flat';
b.CData(2,:) = [0.7 0.3 0.3];
b.CData(3,:) = [0.3 0.7 0.3];
errorbar(x,CR_errorplot,std_CR,'k.')
xticks([1 2 3])
xticklabels({'LR' 'GPR' 'GPRc'})
xlabel('Regression model','FontSize',16)
ylabel('Completion rate')