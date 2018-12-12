%% Representative trajectorie from all models
hfig_traj = figure;
hold on
set(hfig_traj, 'units', 'normalized','renderer','painters','name','Experiment','numbertitle','off','Color','w')

subplot(2,3,1)
hold on
plot(Subject8_GoodLR(:,1),Subject8_GoodLR(:,2),'.b')
%plot(Subject6_BadLR(:,1),Subject6_BadLR(:,2),'.r')
axis([-1 1 -1 1])
set(gca, 'color', [.98,.98,.98],'DataAspectRatio',[1 1 1], 'units', 'normalized', 'xtick', 0, 'ytick', 0, 'LineWidth', 1,'GridColor','k')%,'drawmode', 'fast')
grid on
title('LR','Fontsize',14)

subplot(2,3,2)
hold on
plot(Subject3_GoodGPR(:,1),Subject3_GoodGPR(:,2),'.b')
%plot(Subject5_BadGPR(:,1),Subject5_BadGPR(:,2),'.r')
axis([-1 1 -1 1])
set(gca, 'color', [.98,.98,.98],'DataAspectRatio',[1 1 1], 'units', 'normalized', 'xtick', 0, 'ytick', 0, 'LineWidth', 1,'GridColor','k')%,'drawmode', 'fast')
grid on
title('GPR','Fontsize',14)

subplot(2,3,3)
hold on
plot(Subject9_GoodGPRc(:,1),Subject9_GoodGPRc(:,2),'.b')
%plot(Subject10_BadGPRc(:,1),Subject10_BadGPRc(:,2),'.r')
axis([-1 1 -1 1])
set(gca, 'color', [.98,.98,.98],'DataAspectRatio',[1 1 1], 'units', 'normalized', 'xtick', 0, 'ytick', 0, 'LineWidth', 1,'GridColor','k')%,'drawmode', 'fast')
grid on
title('GPRc','Fontsize',14)

subplot(2,3,4)
hold on
%plot(Subject8_GoodLR(:,1),Subject8_GoodLR(:,2),'.b')
plot(Subject6_BadLR(:,1),Subject6_BadLR(:,2),'.r')
axis([-1 1 -1 1])
set(gca, 'color', [.98,.98,.98],'DataAspectRatio',[1 1 1], 'units', 'normalized', 'xtick', 0, 'ytick', 0, 'LineWidth', 1,'GridColor','k')%,'drawmode', 'fast')
grid on
%title('LR','Fontsize',14)

subplot(2,3,5)
hold on
%plot(Subject3_GoodGPR(:,1),Subject3_GoodGPR(:,2),'.b')
plot(Subject5_BadGPR(:,1),Subject5_BadGPR(:,2),'.r')
axis([-1 1 -1 1])
set(gca, 'color', [.98,.98,.98],'DataAspectRatio',[1 1 1], 'units', 'normalized', 'xtick', 0, 'ytick', 0, 'LineWidth', 1,'GridColor','k')%,'drawmode', 'fast')
grid on
%title('GPR','Fontsize',14)

subplot(2,3,6)
hold on
%plot(Subject9_GoodGPRc(:,1),Subject9_GoodGPRc(:,2),'.b')
plot(Subject10_BadGPRc(:,1),Subject10_BadGPRc(:,2),'.r')
axis([-1 1 -1 1])
set(gca, 'color', [.98,.98,.98],'DataAspectRatio',[1 1 1], 'units', 'normalized', 'xtick', 0, 'ytick', 0, 'LineWidth', 1,'GridColor','k')%,'drawmode', 'fast')
grid on
%title('GPRc','Fontsize',14)