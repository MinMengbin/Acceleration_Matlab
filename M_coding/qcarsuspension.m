
% Physical parameters
mb = 300;    % kg
mw = 60;     % kg
bs = 1000;   % N/m/s
ks = 16000 ; % N/m
kt = 190000; % N/m

% State matrices
A = [ 0 1 0 0; [-ks -bs ks bs]/mb ;0 0 0 1; [ks bs -ks-kt -bs]/mw];
B = [ 0 0; 0 10000/mb ; 0 0 ; [kt -10000]/mw];
C = [1 0 0 0; 1 0 -1 0; A(2,:)];
D = [0 0; 0 0; B(2,:)];

qcar = ss(A,B,C,D);
qcar.StateName = {'body travel (m)';'body vel (m/s)';'wheel travel (m)';'wheel vel (m/s)'};
qcar.InputName = {'r';'fs'};
qcar.OutputName = {'xb';'sd';'bd'};
%zero(qcar('sd','fs'));
%tzero(qcar({'xb','bd'},'fs'));
bodemag(qcar({'bd','sd'},'r'),'b',qcar({'bd','sd'},'fs'),'r',{1 1000});
legend('Road disturbance (r)','Actuator force (fs)','location','SouthWest');
title(['Gain from road dist (r) and actuator force (fs) ' 'to body accel (bd) and suspension travel (sd)']);