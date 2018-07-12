% clear all
% 
% m1 = 2500;
% m2 = 320;
% k1 = 80000;  %k1 = 80000;
% k2 = 500000;
% b1 = 350;
% b2 = 15020; %b2 = 15020;
% 
% kd = 208025;
% kp = 832100;
% ki = 624075;
% 
% 
% 
% nump=[(m1+m2) b2 k2];
% denp=[(m1*m2) (m1*(b1+b2))+(m2*b1) (m1*(k1+k2))+(m2*k1)+(b1*b2) (b1*k2)+(b2*k1) k1*k2];
% G1=tf(nump,denp);
% 
% num1=[-(m1*b2) -(m1*k2) 0 0];
% den1=[(m1*m2) (m1*(b1+b2))+(m2*b1) (m1*(k1+k2))+(m2*k1)+(b1*b2) (b1*k2)+(b2*k1) k1*k2];
% G2=tf(num1,den1);
% 
% numf=num1;
% denf=nump;
% F=tf(numf,denf);
% 
% 
% C = pid(kp,ki,kd);
% 
% sys_cl=F*feedback(G1,C);
% 
% figure;
% t=0:0.05:5;
% step(0.1*sys_cl,t)
% title('Response to a 0.1-m Step under PID Control')

%% x1/w with u  _passive suspension layout

clear all
% parameters setup
parameters_modelv3; % load parameters

% set PID parameters

% kd = 208025;
% kp = 832100;
% ki = 624075;

kd = 0;
kp = 0;
ki = 0;

% laplace transform
s = tf('s');
% force control input x1/U
G10 =  (M2*s^2+b2*s+K2)/((M1*s^2+b1*s+K1)*(M2*s^2+(b1+b2)*s+(K1+K2))-(b1*s+K1)*(b1*s+K1));
% road disturbance input, x1/W
G3 =  (b1*b2*s^2+(b1*K2+b2*K1)*s+K1*K2)/((M1*s^2+b1*s+K1)*(M2*s^2+(b1+b2)*s+(K1+K2))-(b1*s+K1)*(b1*s+K1));

% nump=[(m1+m2) b2 k2];
% denp=[(m1*m2) (m1*(b1+b2))+(m2*b1) (m1*(k1+k2))+(m2*k1)+(b1*b2) (b1*k2)+(b2*k1) k1*k2];
% G1=tf(nump,denp);

[nump,denp] = tfdata(G10);

% num1=[-(m1*b2) -(m1*k2) 0 0];
% den1=[(m1*m2) (m1*(b1+b2))+(m2*b1) (m1*(k1+k2))+(m2*k1)+(b1*b2) (b1*k2)+(b2*k1) k1*k2];
% G2=tf(num1,den1);

[num1,den1] = tfdata(G3);

numf=num1;
denf=nump;
F=tf(numf,denf);


C = pid(kp,ki,kd);

sys_cl=F*feedback(G10,0.1*C);

figure;
t=0:0.05:3;
step(0.1*sys_cl,t)
title('Response to a 0.1-m Step under PID Control')
