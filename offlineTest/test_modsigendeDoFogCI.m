close all

Extension = 1;   % 1 = extension, 2 = flexion, 3 = radial dev., 4 = ulnar dev. 
Flexion = 2;
RadDev = 3;
UlnDev = 4;

N = length(dof);
test = NaN * ones(N,2);
testCI = NaN * ones(N,2);
tl = linspace(0,20,N);

%% X aksen
%% Extension most confidence but higher output on Flexion
for ii = 1:N
   if ci(ii,Extension)<ci(ii,Flexion) && dof(ii,Extension)<dof(ii,Flexion)
       test(ii,1)=dof(ii,Extension);
       test(ii,2)=dof(ii,Flexion);
       testCI(ii,1)=ci(ii,Extension);
       testCI(ii,2)=ci(ii,Flexion);
   end
end
Percentage_llExt = sum(~isnan(test(:,1)))/N * 100;

figure(1)
plot(tl,test,'-or')
hold on; plot(tl,testCI, '--or')

test = NaN * ones(N,2);
testCI = NaN * ones(N,2);
%% Flexion most confidence but higher output on Extension
for ii = 1:N
   if ci(ii,Extension)>ci(ii,Flexion) && dof(ii,Extension)>dof(ii,Flexion)
       test(ii,1)=dof(ii,Extension);
       test(ii,2)=dof(ii,Flexion);
       testCI(ii,1)=ci(ii,Extension);
       testCI(ii,2)=ci(ii,Flexion);
   end
end
Percentage_llFlex = sum(~isnan(test(:,1)))/N * 100;

plot(tl,test,'-or')
plot(tl,testCI, '--or')

test = NaN * ones(N,2);
testCI = NaN * ones(N,2);
%% Extension most confidence and highest output
for ii = 1:N
   if ci(ii,Extension)<ci(ii,Flexion) && dof(ii,Extension)>dof(ii,Flexion)
       test(ii,1)=dof(ii,Extension);
       test(ii,2)=dof(ii,Flexion);
       testCI(ii,1)=ci(ii,Extension);
       testCI(ii,2)=ci(ii,Flexion);
   end
end
Percentage_lhExt = sum(~isnan(test(:,1)))/N * 100;

plot(tl,test,'-og')
plot(tl,testCI, '--og')

test = NaN * ones(N,2);
testCI = NaN * ones(N,2);
%% Flexion most confidence and highest output
for ii = 1:N
   if ci(ii,Extension)>ci(ii,Flexion) && dof(ii,Extension)<dof(ii,Flexion)
       test(ii,1)=dof(ii,Extension);
       test(ii,2)=dof(ii,Flexion);
       testCI(ii,1)=ci(ii,Extension);
       testCI(ii,2)=ci(ii,Flexion);
   end
end
Percentage_lhFlex = sum(~isnan(test(:,1)))/N * 100;

plot(tl,test,'-og')
plot(tl,testCI, '--og')

test = NaN * ones(N,2);
testCI = NaN * ones(N,2);

%% Flexion and Extension have equal confidence or output
% happens due to our confidence treshold "filtration" were some movements is
% zeroed becuase of low confidence
for ii = 1:N
   if ci(ii,Extension)==ci(ii,Flexion) || dof(ii,Extension)==dof(ii,Flexion)
       test(ii,1)=dof(ii,Extension);
       test(ii,2)=dof(ii,Flexion);
       testCI(ii,1)=ci(ii,Extension);
       testCI(ii,2)=ci(ii,Flexion);
   end
end
Percentage_FlexEqualsExt = sum(~isnan(test(:,1)))/N * 100;

plot(tl,test,'-oy')
plot(tl,testCI, '--oy')

test = NaN * ones(N,2);
testCI = NaN * ones(N,2);
%% Y aksen
%% Radial most confidence but higher output on Ulnar
for ii = 1:N
   if ci(ii,RadDev)<ci(ii,UlnDev) && dof(ii,RadDev)<dof(ii,UlnDev)
       test(ii,1)=dof(ii,RadDev);
       test(ii,2)=dof(ii,UlnDev);
       testCI(ii,1)=ci(ii,RadDev);
       testCI(ii,2)=ci(ii,UlnDev);
   end
end
Percentage_llRad = sum(~isnan(test(:,1)))/N * 100;

figure(2)
plot(tl,test,'-or')
hold on; plot(tl,testCI, '--or')

test = NaN * ones(N,2);
testCI = NaN * ones(N,2);
%% Ulnar most confidence but higher output on Radial
for ii = 1:N
   if ci(ii,RadDev)>ci(ii,UlnDev) && dof(ii,RadDev)>dof(ii,UlnDev)
       test(ii,1)=dof(ii,RadDev);
       test(ii,2)=dof(ii,UlnDev);
       testCI(ii,1)=ci(ii,RadDev);
       testCI(ii,2)=ci(ii,UlnDev);
   end
end
Percentage_llUln = sum(~isnan(test(:,1)))/N * 100;

plot(tl,test,'-or')
plot(tl,testCI, '--or')

test = NaN * ones(N,2);
testCI = NaN * ones(N,2);
%% Radial most confidence and highest output
for ii = 1:N
   if ci(ii,RadDev)<ci(ii,UlnDev) && dof(ii,RadDev)>dof(ii,UlnDev)
       test(ii,1)=dof(ii,RadDev);
       test(ii,2)=dof(ii,UlnDev);
       testCI(ii,1)=ci(ii,RadDev);
       testCI(ii,2)=ci(ii,UlnDev);
   end
end
Percentage_lhRad = sum(~isnan(test(:,1)))/N * 100;

plot(tl,test,'-og')
plot(tl,testCI, '--og')

test = NaN * ones(N,2);
testCI = NaN * ones(N,2);
%% Ulnar most confidence and highest output
for ii = 1:N
   if ci(ii,RadDev)>ci(ii,UlnDev) && dof(ii,RadDev)<dof(ii,UlnDev)
       test(ii,1)=dof(ii,RadDev);
       test(ii,2)=dof(ii,UlnDev);
       testCI(ii,1)=ci(ii,RadDev);
       testCI(ii,2)=ci(ii,UlnDev);
   end
end
Percentage_lhUln = sum(~isnan(test(:,1)))/N * 100;

plot(tl,test,'-og')
plot(tl,testCI, '--og')

test = NaN * ones(N,2);
testCI = NaN * ones(N,2);

%% Radial and Ulna have equal confidence or output
% happens due to our confidence treshold "filtration" were some movements is
% zeroed becuase of low confidence
for ii = 1:N
   if ci(ii,RadDev)==ci(ii,UlnDev) || dof(ii,RadDev)==dof(ii,UlnDev)
       test(ii,1)=dof(ii,RadDev);
       test(ii,2)=dof(ii,UlnDev);
       testCI(ii,1)=ci(ii,RadDev);
       testCI(ii,2)=ci(ii,UlnDev);
   end
end
Percentage_RadEqualsUln = sum(~isnan(test(:,1)))/N * 100;

plot(tl,test,'-oy')
plot(tl,testCI, '--oy')

test = NaN * ones(N,2);
testCI = NaN * ones(N,2);
