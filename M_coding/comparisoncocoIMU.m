load('C:\Users\mm395\Documents\MATLAB\REC2098.mat'); % load acceleration data
load('C:\Users\mm395\Documents\MATLAB\Untitled.mat'); % load acceleration data
plot((Untitled(:,1)-Untitled(1,1))/1e9,Untitled(:,2),'DisplayName','Untitled');
hold on 
plot(REC2098(:,1)+2,REC2098(:,2),'-');
hold off;