 %%%%%  load data %%%%%%
REC0826ch2 = load('C:\Users\mm395\Documents\MATLAB\REC0826ch2.mat'); % load acceleration data

%%%%%  FFT analysis of raw data %%%%%%
figure;
subplot(2,1,1);
plot(REC0826ch2.REC0826ch2(:,1),REC0826ch2.REC0826ch2(:,2));
hold on
title('Acceleration data');
fft_accel = fft(REC0826ch2.REC0826ch2(:,2)); % Compute the FFT of the acceleration data
L_accel = length(REC0826ch2.REC0826ch2(:,2)); % Calculate the number of the data
P_2 = abs(fft_accel/L_accel); % Compute the two-sided spectrum P2
P_1 = P_2(1:L_accel/2+1); % Compute the sing-sided spectrum P1
P_1(2:end-1) = 2*P_1(2:end-1);
Fs = 3200; %1/(137/L_accel);%
f_accel = Fs*(0:(L_accel/2))/L_accel;
subplot(2,1,2);
plot(f_accel,P_1); 
title('Single-Sided Amplitude Spectrum of Acceleration Data');
xlabel('f (Hz)');
ylabel('|P1(f)|'); 

%%%%%  filter data %%%%%%
% %apply low pass filter
% fc = 100;  %cut off frequency
% fs = Fs;
% [b,a] = butter(6,fc/(fs/2));
% freqz(b,a);
% 
% %apply high pass filter
% [z,p,k] = butter(9,300/500,'high');
% sos = zp2sos(z,p,k);
% fvtool(sos,'Analysis','freq')

% %apply band pass filter
d = fdesign.bandpass('N,F3dB1,F3dB2',10,0.5,50,3200);
Hd = design(d,'butter');

% Frequency analysis of Filtered acceleration data
% filtered_accel = filter(b, a, REC0826ch2(:,2) ); % low pass

filtered_accel = filter(Hd, REC0826ch2.REC0826ch2(:,2) ); % band pass

figure;plot(REC0826ch2.REC0826ch2(:,1),filtered_accel);
title(' Filtered Acceleration data');
fft_accel_f = fft(filtered_accel); % Compute the FFT of the acceleration data
P_2_f = abs(fft_accel_f/L_accel); % Compute the two-sided spectrum P2
P_1_f = P_2_f(1:L_accel/2+1); % Compute the sing-sided spectrum P1
P_1_f(2:end-1) = 2*P_1_f(2:end-1);
Fs = 3200; %1/(137/L_accel);%
f_accel_f = Fs*(0:(L_accel/2))/L_accel;
figure;plot(f_accel_f,P_1_f); 
title('Single-Sided Amplitude Spectrum of Filtered Acceleration Data');
xlabel('f (Hz)');
ylabel('|P1(f)|'); 

%%%%%  caculate the accumulative velocity (integration if acceleration data) %%%%%%
V_acce = trapz(filtered_accel);
AC_V_acce = cumtrapz(filtered_accel)*(1/Fs);
figure; plot(REC0826ch2.REC0826ch2(:,1),AC_V_acce);
title('Speed');
xlabel('Time (s)');
ylabel('Velocity (m/s)');

%%%%%  FFT analysis of velocity data %%%%%%
% fft_velo = fft(AC_V_acce); % Compute the FFT of the acceleration data
% L_velo = length(AC_V_acce); % Calculate the number of the data
% P_2 = abs(fft_velo/L_velo); % Compute the two-sided spectrum P2
% P_1 = P_2(1:L_velo/2+1); % Compute the sing-sided spectrum P1
% P_1(2:end-1) = 2*P_1(2:end-1);
% f_velo = Fs*(0:(L_velo/2))/L_velo;
% figure;plot(f_velo,P_1); 
% title('Single-Sided Amplitude Spectrum of velocity Data');
% xlabel('f (Hz)');
% ylabel('|P1(f)|'); 


% %apply low pass filter
% filtered_V_acce = filter(b, a, V_acce);
% 
% 
%caculate accumulative displacement
% D_acce = trapz(filtered_V_acce); % filtered velocity
% AD_V_acce = cumtrapz(filtered_V_acce)*(1/Fs); % filtered velocity
D_acce = trapz(AC_V_acce);
AD_V_acce = cumtrapz(AC_V_acce)*(1/Fs);
% T = table(time',AC_V_acce','VariableNames',{'Time','CumulativeDistance'});
figure; plot(REC0826ch2.REC0826ch2(:,1),AD_V_acce);
title('Displacement');
xlabel('Time (s)');
ylabel('Displacement (m)');


% plot the result
% figure; plot(REC0826ch2.REC0826ch2(:,1),REC0826ch2.REC0826ch2(:,2), 'r');
% title('Comparision between filted acceleration and raw accleration');
% xlabel('Time (s)');
% ylabel('Acceleration (m/s^2)');
% hold on
% plot(REC0826ch2.REC0826ch2(:,1),filtered_accel,'k');


%creat timeseries of displacement data

% REC0826ch2_dis = timeseries(AD_V_acce,REC0826ch2.REC0826ch2(:,1));
