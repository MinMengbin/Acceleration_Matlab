 %%%%%  load data Engin off data %%%%%%
REC0804ch1 = load('C:\Users\mm395\Documents\MATLAB\REC0804ch1.mat'); % load acceleration data

%%  FFT analysis of raw data %%%%%%
figure;
subplot(2,1,1);
plot(REC0804ch1.REC0804ch1(:,1),REC0804ch1.REC0804ch1(:,2));
title('Acceleration data');
ylabel('Acceleration m/s/s');
xlabel('Time /second'); 
hold on
fft_accel = fft(REC0804ch1.REC0804ch1(:,2)); % Compute the FFT of the acceleration data
L_accel = length(REC0804ch1.REC0804ch1(:,2)); % Calculate the number of the data
P_2 = abs(fft_accel/L_accel); % Compute the two-sided spectrum P2
P_1 = P_2(1:L_accel/2+1); % Compute the sing-sided spectrum P1
P_1(2:end-1) = 2*P_1(2:end-1);
Fs = 800; %1/(137/L_accel);%
f_accel = Fs*(0:(L_accel/2))/L_accel;
subplot(2,1,2);
plot(f_accel,P_1); 
title('Single-Sided Amplitude Spectrum of Acceleration Data');
xlabel('f (Hz)');
ylabel('|P1(f)|'); 
hold on 
% %% apply band pass filter
% d = fdesign.bandpass('N,F3dB1,F3dB2',10,0.001,1000,3200);
% Hd = design(d,'butter');
% 
% %% % Analysis of Filtered acceleration data
% 
% filtered_accel = filter(Hd, REC0804ch1.REC0804ch1(:,2) );
% subplot(2,1,1);
% plot(REC0804ch1.REC0804ch1(:,1),filtered_accel);
% fft_accel_f = fft(filtered_accel); % Compute the FFT of the acceleration data
% P_2_f = abs(fft_accel_f/L_accel); % Compute the two-sided spectrum P2
% P_1_f = P_2_f(1:L_accel/2+1); % Compute the sing-sided spectrum P1
% P_1_f(2:end-1) = 2*P_1_f(2:end-1);
% Fs = 3200; %1/(137/L_accel);%
% f_accel_f = Fs*(0:(L_accel/2))/L_accel;
% subplot(2,1,2);
% plot(f_accel_f,P_1_f); 
% hold off