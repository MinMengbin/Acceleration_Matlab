openfig('matlab60.fig');
load('C:\Users\mm395\Documents\MATLAB\REC2100.mat'); % load acceleration data
hold on
plot(REC2100(:,1)+0.26,REC2100(:,2),'r');
hold off