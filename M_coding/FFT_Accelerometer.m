
% %test coding for FFT
% t = 0:.001:1.499;
% % for i = 1 : 1500
% %  x(i) = 2;
% % end
% %x = 1*sin(2*pi*120*t);
% figure;plot(t,x);
% title('Acceleration data');
% fft_accel = fft(x); % Compute the FFT of the acceleration data
% L_accel = length(x); % Calculate the number of the data
% P_2 = abs(fft_accel/L_accel); % Compute the two-sided spectrum P2
% P_1 = P_2(1:L_accel/2+1); % Compute the sing-sided spectrum P1
% P_1(2:end-1) = 2*P_1(2:end-1);
% Fs = 1000; %1/(137/L_accel);%
% f_accel = Fs*(0:(L_accel/2))/L_accel;
% figure;plot(f_accel,P_1); 
% title('Single-Sided Amplitude Spectrum of X(t)');
% xlabel('f (Hz)');
% ylabel('|P1(f)|'); 



load('C:\Users\mm395\Documents\MATLAB\REC0826ch2.mat'); % load acceleration data
figure;plot(REC0826ch2(:,1),REC0826ch2(:,2));
title('Acceleration data');
fft_accel = fft(REC0826ch2(:,2)); % Compute the FFT of the acceleration data
L_accel = length(REC0826ch2(:,2)); % Calculate the number of the data
P_2 = abs(fft_accel/L_accel); % Compute the two-sided spectrum P2
P_1 = P_2(1:L_accel/2+1); % Compute the sing-sided spectrum P1
P_1(2:end-1) = 2*P_1(2:end-1);
Fs = 3200; %1/(137/L_accel);%
f_accel = Fs*(0:(L_accel/2))/L_accel;
figure;plot(f_accel,P_1); 
title('Single-Sided Amplitude Spectrum of X(t)');
xlabel('f (Hz)');
ylabel('|P1(f)|'); 


% %%% coding used to test the cumtrapz function
% t = 0:0.1:1.9;
% for i = 1 : 20
%  x(i) = 2;
% end
% figure;plot(t,x);
% title('Acceleration data');
% 
% AC_V_acce = cumtrapz(x)*0.1;
% figure;plot(t,AC_V_acce);
% title('velocity');
% %%% coding used to test the cumtrapz function



%apply filter of the data

%caculate the accumulative velocity
time = 91.2; % Data recording time
V_acce = trapz(REC0826ch2(:,2));
AC_V_acce = cumtrapz(REC0826ch2(:,2))*(1/Fs);
%T = table(time',AC_V_acce','VariableNames',{'Time','CumulativeDistance'});
figure; plot(REC0826ch2(:,1),AC_V_acce);
title('Speed');
xlabel('Time (s)');
ylabel('Velocity (m/s^2)');

% %acceleration integration to velocity 
% V_acce = zeros(L_accel);
% 
% %V_acce = trapz(x);
% for i = 1:L_accel-1
%   V_acce(i+1) = ((x(i+1) + x(i))/2)* (t(i+1) - t(i));
% end
% %T = table(time',AC_V_acce','VariableNames',{'Time','CumulativeDistance'});
% figure; plot(t,V_acce);
% title('Speed');
% xlabel('Time (s)');
% ylabel('Velocity (m/s)');



% % test coding for acceleration integration to velocity 
% V_acce = zeros(L_accel,1);
% for i = 1:L_accel-1
%   V_acce(i+1) = ((REC0826ch2(i+1,2) + REC0826ch2(i,2))/2) * (REC0826ch2(i+1,1) - REC0826ch2(i,1));
% end
% %T = table(time',AC_V_acce','VariableNames',{'Time','CumulativeDistance'});
% figure; plot(REC0826ch2(:,1),V_acce);
% title('Speed');
% xlabel('Time (s)');
% ylabel('Velocity (m/s)');




%caculate accumulative displacement
D_acce = trapz(AC_V_acce);
AD_V_acce = cumtrapz(AC_V_acce)*(1/Fs);
%T = table(time',AC_V_acce','VariableNames',{'Time','CumulativeDistance'});
figure; plot(REC0826ch2(:,1),AD_V_acce);
title('Displacement');
xlabel('Time (s)');
ylabel('Displacement (m)');


% % test coding for velocity integration to displacement
% D_acce = zeros(L_accel);
% 
% %V_acce = trapz(x);
% for i = 1:L_accel-1
%   D_acce(i+1) = ((V_acce(i+1) + V_acce(i))/2)* (t(i+1) - t(i));
% end
% %T = table(time',AC_V_acce','VariableNames',{'Time','CumulativeDistance'});
% figure; plot(t,D_acce);
% title('Displacement');
% xlabel('Time (s)');
% ylabel('Velocity (m)');

% % test coding for acceleration integration to velocity 
% D_acce = zeros(L_accel,1);
% 
% for i = 1:L_accel-1
%   D_acce(i+1) = ((V_acce(i+1) + V_acce(i))/2)* (REC0826ch2(i+1,1) - REC0826ch2(i,1));
% end
% figure; plot(REC0826ch2(:,1),D_acce);
% title('Displacement');
% xlabel('Time (s)');
% ylabel('Displacement (m)');
