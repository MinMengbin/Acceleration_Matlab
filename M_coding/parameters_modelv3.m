%%% this coding is used to change the parameters' value for the
clear all;
%% Suspension_tractor_V2.slx model
kt = 1.074383224580817e+05; %  the stiffness of the tyre is 104420 n/m (from experiments)
K2 = kt; % variable name for passive_suspension_ss.m
k2 = kt; % variable name for passive_suspension_ss.m
ks = 100000; %  the constant of spring is 2000n/m
K1 = ks; % variable name for passive_suspension_ss.m
k1 = ks; % variable name for passive_suspension_ss.m

dc = 1000 ; %  the coefficient of the damper is 6000 n*s/m
b1 = dc; % variable name for passive_suspension_ss.m
b2 = 0; % the coefficient of tyre dampness

m = 150/4; % frame mass is around 150 kg, but this is a quarter system model
M1 = m; % variable name for passive_suspension_ss.m
m1 = m; % variable name for passive_suspension_ss.m
mu = 20/2; % unprund frame mass is 30 kg
M2 = mu; % variable name for passive_suspension_ss.m
m2 = m; % variable name for passive_suspension_ss.m

