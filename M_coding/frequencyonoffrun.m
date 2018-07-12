 %%%%%  load data Engin off data %%%%%%
REC0804ch1 = load('C:\Users\mm395\Documents\MATLAB\REC0804ch1.mat'); % load acceleration data

%%  FFT analysis of raw data %%%%%%
figure;
subplot(2,3,1);
plot(REC0804ch1.REC0804ch1(:,1),REC0804ch1.REC0804ch1(:,2));
title('Acceleration data when engin is off');
ylabel('Acceleration m/s/s');
xlabel('Time /second'); 
hold on
fft_accel = fft(REC0804ch1.REC0804ch1(:,2)); % Compute the FFT of the acceleration data
L_accel = length(REC0804ch1.REC0804ch1(:,2)); % Calculate the number of the data
P_2 = abs(fft_accel/L_accel); % Compute the two-sided spectrum P2
P_1 = P_2(1:L_accel/2+1); % Compute the sing-sided spectrum P1
P_1(2:end-1) = 2*P_1(2:end-1);
Fs = 1/0.0001;
f_accel = Fs*(0:(L_accel/2))/L_accel;
subplot(2,3,4);
plot(f_accel,P_1); 
title('Single-Sided Amplitude Spectrum of Acceleration Data');
xlabel('f (Hz)');
ylabel('|P1(f)|'); 
hold on 
%% apply band pass filter
d = fdesign.bandpass('N,F3dB1,F3dB2',10,0.001,1000,3200);
Hd = design(d,'butter');

%% % Analysis of Filtered acceleration data

filtered_accel = filter(Hd, REC0804ch1.REC0804ch1(:,2) );
subplot(2,3,1);
plot(REC0804ch1.REC0804ch1(:,1),filtered_accel);
fft_accel_f = fft(filtered_accel); % Compute the FFT of the acceleration data
P_2_f = abs(fft_accel_f/L_accel); % Compute the two-sided spectrum P2
P_1_f = P_2_f(1:L_accel/2+1); % Compute the sing-sided spectrum P1
P_1_f(2:end-1) = 2*P_1_f(2:end-1);
Fs = 1/0.0001; 
f_accel_f = Fs*(0:(L_accel/2))/L_accel;
subplot(2,3,4);
plot(f_accel_f,P_1_f); 
hold off



 %%%%%  load data of Engine on %%%%%%
REC0810ch2 = load('C:\Users\mm395\Documents\MATLAB\REC0810ch2.mat'); % load acceleration data

%%  FFT analysis of raw data %%%%%%

subplot(2,3,2);
plot(REC0810ch2.REC0810ch2(:,1),REC0810ch2.REC0810ch2(:,2));
title('Acceleration data when engin is on');
ylabel('Acceleration m/s/s');
xlabel('Time /second'); 
hold on
fft_accel = fft(REC0810ch2.REC0810ch2(:,2)); % Compute the FFT of the acceleration data
L_accel = length(REC0810ch2.REC0810ch2(:,2)); % Calculate the number of the data
P_2 = abs(fft_accel/L_accel); % Compute the two-sided spectrum P2
P_1 = P_2(1:L_accel/2+1); % Compute the sing-sided spectrum P1
P_1(2:end-1) = 2*P_1(2:end-1);
Fs = 3200; %1/(137/L_accel);%
f_accel = Fs*(0:(L_accel/2))/L_accel;
subplot(2,3,5);
plot(f_accel,P_1); 
title('Single-Sided Amplitude Spectrum of Acceleration Data');
xlabel('f (Hz)');
ylabel('|P1(f)|'); 
hold on 

%% apply band pass filter
d = fdesign.bandpass('N,F3dB1,F3dB2',10,0.001,1000,3200);
Hd = design(d,'butter');

%% % Analysis of Filtered acceleration data

filtered_accel = filter(Hd, REC0810ch2.REC0810ch2(:,2) );
subplot(2,3,2);
plot(REC0810ch2.REC0810ch2(:,1),filtered_accel);
fft_accel_f = fft(filtered_accel); % Compute the FFT of the acceleration data
P_2_f = abs(fft_accel_f/L_accel); % Compute the two-sided spectrum P2
P_1_f = P_2_f(1:L_accel/2+1); % Compute the sing-sided spectrum P1
P_1_f(2:end-1) = 2*P_1_f(2:end-1);
Fs = 3200; %1/(137/L_accel);%
f_accel_f = Fs*(0:(L_accel/2))/L_accel;
subplot(2,3,5);
plot(f_accel_f,P_1_f);  
hold off

 %%%%%  load data of Engine on %%%%%%
REC0826ch2 = load('C:\Users\mm395\Documents\MATLAB\REC0826ch2.mat'); % load acceleration data

%%  FFT analysis of raw data %%%%%%

subplot(2,3,3);

plot(REC0826ch2.REC0826ch2(:,1),REC0826ch2.REC0826ch2(:,2));
title('Acceleration data when it is running');
ylabel('Acceleration m/s/s');
xlabel('Time /second'); 
hold on
fft_accel = fft(REC0826ch2.REC0826ch2(:,2)); % Compute the FFT of the acceleration data
L_accel = length(REC0826ch2.REC0826ch2(:,2)); % Calculate the number of the data
P_2 = abs(fft_accel/L_accel); % Compute the two-sided spectrum P2
P_1 = P_2(1:L_accel/2+1); % Compute the sing-sided spectrum P1
P_1(2:end-1) = 2*P_1(2:end-1);
Fs = 3200; %1/(137/L_accel);%
f_accel = Fs*(0:(L_accel/2))/L_accel;
subplot(2,3,6);
plot(f_accel,P_1); 
title('Single-Sided Amplitude Spectrum of Acceleration Data');
xlabel('f (Hz)');
ylabel('|P1(f)|'); 

%% apply band pass filter
d = fdesign.bandpass('N,F3dB1,F3dB2',10,0.001,1000,3200);
Hd = design(d,'butter');

%% % Analysis of Filtered acceleration data
hold on 
filtered_accel = filter(Hd, REC0826ch2.REC0826ch2(:,2) );
subplot(2,3,3);
plot(REC0826ch2.REC0826ch2(:,1),filtered_accel);

fft_accel_f = fft(filtered_accel); % Compute the FFT of the acceleration data
P_2_f = abs(fft_accel_f/L_accel); % Compute the two-sided spectrum P2
P_1_f = P_2_f(1:L_accel/2+1); % Compute the sing-sided spectrum P1
P_1_f(2:end-1) = 2*P_1_f(2:end-1);
Fs = 3200; %1/(137/L_accel);%
f_accel_f = Fs*(0:(L_accel/2))/L_accel;
subplot(2,3,6);
plot(f_accel_f,P_1_f); 

hold off