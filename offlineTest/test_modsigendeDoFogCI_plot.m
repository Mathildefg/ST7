%% Bar plot of contrary movement percentage

S = std(data);
m = mean(data);
Y = [m(2),m(3),m(6),m(7)
     m(4),m(5),m(8),m(9)];

Y_errorplot = [m(2);m(3);m(6);m(7);m(4);m(5);m(8);m(9)];
x = [0.725; 0.91; 1.09; 1.275;
    1.725; 1.91; 2.09; 2.275]; 
 
std_Y =[S(2);S(3);S(6);S(7);S(4);S(5);S(8);S(9)];
 
figure
set(gca,'FontSize',18) % Creates an axes and sets its FontSize to 18
hold on
bar(Y)
errorbar(x,Y_errorplot,std_Y,'.')
xticks([1 2])
xticklabels({'Extension vs. Flexion' 'Radial vs. Ulnar'})
xlabel('Movements axis')
yticks([10 20 30 40 50 60 70 80 90 100])
yticklabels({'10%' '20%' '30%' '40%'})
ylabel('Percentage of signal')
title('Contradicting Predictions')
legend('Low CI and high output','High CI and low output','Low CI and low output','High CI and high output','STD')
