load('C:\Users\mm395\Documents\MATLAB\REC0826ch1.mat'); % load acceleration data
figure;plot(REC0826ch1(:,1),REC0826ch1(:,2));
title('Displacement data');
xlabel('time');
ylabel('m'); 
fft_accel = fft(REC0826ch1(:,2)); % Compute the FFT of the acceleration data
L_accel = length(REC0826ch1(:,2)); % Calculate the number of the data
P_2 = abs(fft_accel/L_accel); % Compute the two-sided spectrum P2
P_1 = P_2(1:L_accel/2+1); % Compute the sing-sided spectrum P1
P_1(2:end-1) = 2*P_1(2:end-1);
Fs = 3200; %1/(137/L_accel);%
f_accel = Fs*(0:(L_accel/2))/L_accel;
figure;plot(f_accel,P_1); 
title('Single-Sided Amplitude Spectrum of Displacement Data');
xlabel('f (Hz)');
ylabel('|P1(f)|'); 


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

% % %apply band pass filter
% d = fdesign.bandpass('N,F3dB1,F3dB2',10,1,100,3200);
% Hd = design(d,'butter');
% 
% % Frequency analysis of Filtered acceleration data
% %filtered_accel = filter(b, a, REC0826ch1(:,2) );
% % 
% filtered_accel = filter(Hd, REC0826ch1(:,2) );
% 
% figure;plot(REC0826ch1(:,1),filtered_accel);
% title(' Filtered Acceleration data');
% fft_accel_f = fft(filtered_accel); % Compute the FFT of the acceleration data
% P_2_f = abs(fft_accel_f/L_accel); % Compute the two-sided spectrum P2
% P_1_f = P_2_f(1:L_accel/2+1); % Compute the sing-sided spectrum P1
% P_1_f(2:end-1) = 2*P_1_f(2:end-1);
% Fs = 3200; %1/(137/L_accel);%
% f_accel_f = Fs*(0:(L_accel/2))/L_accel;
% figure;plot(f_accel_f,P_1_f); 
% title('Single-Sided Amplitude Spectrum of Filtered Acceleration Data');

