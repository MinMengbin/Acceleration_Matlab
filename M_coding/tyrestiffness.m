clear;
% 
Fs = 25600;
% 
% figure;
% % %% load free drop data
% % Height = 0.07; % 0.07 m
% % M_w = 5.9; % mass of wheel 5.9 kg
 REC2119ch1 = load('C:\Users\mm395\Documents\MATLAB\REC2119ch1.mat'); % load acceleration data
% subplot(3,1,1);
% plot(REC2117ch1.REC2117ch1(3.7*Fs:5.1*Fs,1),REC2117ch1.REC2117ch1(3.7*Fs:5.1*Fs,2));
% title('Raw Acceleration data');
% ylabel('Acceleration m/s/s');
% xlabel('Time /s');
% 
% % % %apply band pass filter
% % % d = fdesign.bandpass('N,F3dB1,F3dB2',10,0.1,1600,3200);
% % % Hd = design(d,'butter');
% % 
% % % filtered_accel = filter(Hd, REC2117ch1.REC2117ch1(3.7*Fs:5.1*Fs,2) ); % band pass
% % % subplot(4,1,2);
% % % plot(REC2117ch1.REC2117ch1(3.7*Fs:5.1*Fs,1),filtered_accel);
% % 
% % 
% filtered_accel = REC2117ch1.REC2117ch1(3.7*Fs:5.1*Fs,2); % no filter applied
% 
% V_acce = trapz(filtered_accel);
% AC_V_acce = cumtrapz(filtered_accel)*(1/Fs);
% subplot(3,1,2); 
% plot(REC2117ch1.REC2117ch1(3.7*Fs:5.1*Fs,1),AC_V_acce);
% title('Speed');
% xlabel('Time (s)');
% ylabel('Velocity (m/s)');
% 
% % AC_V_acce  = filter(Hd, AC_V_acce); % band pass
% 
% 
% D_acce = trapz(AC_V_acce);
% AD_V_acce = cumtrapz(AC_V_acce)*(1/Fs);
% subplot(3,1,3); 
% plot(REC2117ch1.REC2117ch1(3.7*Fs:5.1*Fs,1),AD_V_acce);
% title('Displacement');
% xlabel('Time (s)');
% ylabel('Displacement (m)');


%  load tyre data, FFT analysis, natural frequency %%%%%%
% REC2111ch1 = load('C:\Users\mm395\Documents\MATLAB\REC2111ch1.mat'); % load acceleration data

%% draw raw data and plot the decay
figure;
subplot(2,1,1);
% plot(REC2111ch1.REC2111ch1(3*800:3.2*800,1),REC2111ch1.REC2111ch1(3*800:3.2*800,2));
plot(REC2119ch1.REC2119ch1(4.3*Fs:4.6*Fs,1),REC2119ch1.REC2119ch1(4.3*Fs:4.6*Fs,2));
title('Acceleration data');
ylabel('Acceleration m/s/s');
xlabel('Time /second'); 
hold on
x_posi = [4.338,4.382,4.425,4.468,4.512,4.554,4.597]; % get the points based on the plot
y_posi = [2.703,1.915,1.360,1.022,0.7636,0.5853,0.4587]; %get the points baded on the plot
plot(x_posi,y_posi); % draw positive part of the decay
hold on
x_neg = [4.315, 4.359, 4.403, 4.447, 4.49, 4.533, 4.576];
y_neg = [-3.17, -2.245, -1.639, -1.135, -0.8776, -0.6276, -0.4151];
plot(x_neg, y_neg); % draw negative part of the decay
hold on

%%  FFT analysis of raw data %%%%%%
% fft_accel = fft(REC2119ch1.REC2119ch1(4.3*Fs:4.6*Fs,2)); % Compute the FFT of the acceleration data
% L_accel = length(REC2119ch1.REC2119ch1(4.3*Fs:4.6*Fs,2)); % Calculate the number of the data
% P_2 = abs(fft_accel/L_accel); % Compute the two-sided spectrum P2
% P_1 = P_2(1:L_accel/2+1); % Compute the sing-sided spectrum P1
% P_1(2:end-1) = 2*P_1(2:end-1);
% Fs = 25600; %1/(137/L_accel);%
% f_accel = Fs*(0:(L_accel/2))/L_accel;
% subplot(2,1,2);
% plot(f_accel,P_1); 
% title('Single-Sided Amplitude Spectrum of Acceleration Data');
% xlabel('f (Hz)');
% ylabel('|P1(f)|'); 
% hold on 

%% Power Spectral Density by using the complex conjugate (CONJ) to look the frequency components of signal
L_accel = length(REC2119ch1.REC2119ch1(4.3*Fs:4.6*Fs,2));
fft_accel = fft(REC2119ch1.REC2119ch1(4.3*Fs:4.6*Fs,2),L_accel); % Compute the FFT of the acceleration data
PSD_fft_accel = fft_accel.*conj(fft_accel)/L_accel;
f = Fs/L_accel*(0:L_accel/2-1);
subplot(2,1,2);
plot(f,PSD_fft_accel(1:L_accel/2));
title('Single-Sided PSD Amplitude of Acceleration Data');
xlabel('f (Hz)');
ylabel('|PSD(f)|'); 

%% find damped frequency based on FFT analysis
fd = 23.33; % Hz

%% find damping ratio based on log decrement
x_ln = (1:1:length(x_posi));
y_ln = (1:1:length(x_posi));
for i = 1:length(y_posi)
  y_ln(i) = log(y_posi(1)/y_posi(i));
end
figure
plot (x_ln, y_ln);
title('log decrement');
xlabel('n');
ylabel('ln(A(0)/A(n))'); 

%%
% run the program first and plot the result of the log decrement
% then turn to the figure by using matlab funciton of basic fitting
% apply linear fit and get the slope of the log decrement
% then get the value of the slop which is 
p1_slope = 0.29532;
% equation: slop value = 2*pi*zeta/sqrt(1-zeta^2)
% then get the damping ratio zeta
zeta_ratio = p1_slope/sqrt(4*pi^2+p1_slope^2);

%% calculate natural frequency based on damped frequency and damping ratio
% equation: natural frequency*sqrt(1-zeta^2) = damped frequency
% then get the natural frequency
fn = fd/sqrt(1-zeta_ratio^2); % Hz

% get tyre stiffness kt and dampness d_t based on the value
m = 5; % the mass of the hub
kt = (fn^2)*4*(pi^2)*m; % tyre stiffness 1.04420e+05 n/m
d_tyre = zeta_ratio*2*sqrt(kt*m); % actural damping = damping ratio * critical damping
 
% %% Tyre model x/W  Frequency Response, Time response
% %% tyre model parameters
% s = tf('s');
% G = kt/(m*s^2+2*kt); %% tyre model; 
% figure(1);
%  w = logspace(-1,3,10000);
% frequency response of No suspension system
%     G = kt/(m*s^2+2*kt); %% tyre model;
%     [A, B, C, D] = ssdata(G);
%     e = eig(A);
%     fs = abs(e);
%     [mag,phase] = bode(G,w);
%     subplot(3,1,1);
%     loglog(w, squeeze(mag));
%     title('Bode Diagram of tyre model system');
%     xlabel('Normalized Frequency rad/s'); % x-axis label
%     ylabel('Magnitutide dB'); % y-axis label
%     text(fs(1,1),20,['natural frequency is: ' num2str(fs(1,1)/(2*pi),'%2.3f') 'Hz'] ,'Color','black','FontSize',14);   
%     subplot(3,1,2);
%     semilogx(w, squeeze(phase));
%     xlabel('Normalized Frequency rad/s'); % x-axis label
%     ylabel('Phase degree'); % y-axis label
%     subplot(3,1,3);
%     step(G);
% % %% Tyre model x/W  Frequency Response, Time response

% hold on 
% %% apply band pass filter
% d = fdesign.bandpass('N,F3dB1,F3dB2',10,0.001,1000,3200);
% Hd = design(d,'butter');
% 
% %% % Analysis of Filtered acceleration data
% 
% filtered_accel = filter(Hd, REC2111ch1.REC2111ch1(:,2) );
% subplot(2,1,1);
% plot(REC2111ch1.REC2111ch1(:,1),filtered_accel);
% fft_accel_f = fft(filtered_accel); % Compute the FFT of the acceleration data
% P_2_f = abs(fft_accel_f/L_accel); % Compute the two-sided spectrum P2
% P_1_f = P_2_f(1:L_accel/2+1); % Compute the sing-sided spectrum P1
% P_1_f(2:end-1) = 2*P_1_f(2:end-1);
% Fs = 3200; %1/(137/L_accel);%
% f_accel_f = Fs*(0:(L_accel/2))/L_accel;
% subplot(2,1,2);
% plot(f_accel_f,P_1_f); 
% hold off