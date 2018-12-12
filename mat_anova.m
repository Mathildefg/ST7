%% 1-way RM ANOVA
%completionRate= [];

LR= completionRate(:,1);
GPR= completionRate(:,2);
GPRc= completionRate(:,3);

data=[LR;GPR;GPRc];
cat = {'LR';'LR';'LR';'LR';'LR';'LR';'LR';'LR';'LR';'LR';'GPR';'GPR';'GPR';'GPR';'GPR';'GPR';'GPR';'GPR';'GPR';'GPR';'GPRc';'GPRc';'GPRc';'GPRc';'GPRc';'GPRc';'GPRc';'GPRc';'GPRc';'GPRc'};

t = table(cat,data, 'VariableNames',{'cat','completionRate'});
CR = table([1]','VariableNames',{'Measurements'});

rm = fitrm(t,'completionRate~cat','WithinDesign',CR);
ranovatbl = ranova(rm)

d = mean(completionRate,1)
