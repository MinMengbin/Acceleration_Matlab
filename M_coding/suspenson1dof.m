%% this coding is used for simulation of a passive suspension system one degree of freedom
% reference: http://ctms.engin.umich.edu/CTMS/index.php?example=Suspension&section=SystemModeling
clear all;

%% Set up
parameters_modelv3; % load parameters

%% transfer function
s = tf('s');
    
% X1/W Comparison Between No Suspension, Passive, Skyhook, Groundhook 
% % both Frequency Response and Time Domain)

b1 = 2*1*sqrt(kt*m);  

G5 = 1/(1+M1*s^2/K2);  % Model for No suspension
%%% Model for Passive suspension
G6 = (k2+b1*s)/(m1*s^2+b1*s+k2); % model for passive

G7 = k2/(m1*s^2+b1*s+k2); %  Model for Skyhook 
figure;
% w = logspace(-1,2,100);
w = (0.01:0.01:50*2*pi);
 
% options = bodeoptions;
% options.FreqUnits = 'Hz'; % or 'rad/second', 'rpm', etc.
% frequency response of No suspension system



%bode diagram
bode(G5,w);
hold on
bode(G6,w);
hold on
bode(G7,w);
hold on
legend('No Suspesion','Passive Suspension','Ideal Skyhook');
% 
% for i = 0.1 : 0.2: 1.3
% 
% %bode diagram
% b1 = 2*i*sqrt(kt*m); 
% G7 = k2/(m1*s^2+b1*s+k2); %  Model for Skyhook 
% bode(G7,w);
% hold on
% 
% end 

hold off
%bode diagram
fs = (0.5/pi)*sqrt(kt/m);


