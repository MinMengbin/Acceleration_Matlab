
load('C:\Users\mm395\Documents\MATLAB\REC2097_stable_ch2.mat'); % load integration data
figure;plot(REC2097stable(:,1),REC2097stable(:,2));
title('Acceleration data');
fft_accel = fft(REC2097stable(:,2)); % Compute the FFT of the acceleration data
L_accel = length(REC2097stable(:,2)); % Calculate the number of the data
P_2 = abs(fft_accel/L_accel); % Compute the two-sided spectrum P2
P_1 = P_2(1:L_accel/2+1); % Compute the sing-sided spectrum P1
P_1(2:end-1) = 2*P_1(2:end-1);
Fs = 12800; %1/(137/L_accel);%
f_accel = Fs*(0:(L_accel/2))/L_accel;
figure;plot(f_accel,P_1); 
title('Single-Sided Amplitude Spectrum of Acceleration Data');
xlabel('f (Hz)');
ylabel('|P1(f)|'); 
