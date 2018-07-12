% figure; plot(REC0826ch1(:,1),AD_V_acce, 'r');
% hold on
% plot(REC0826ch1(:,1),AD_100_V_acce,'b');
% %hold on
% plot(REC0826ch1(:,1),AD_1k_V_acce,'y');
% plot(REC0826ch1(:,1),REC0826ch1(:,2),'k');
% 
% figure; plot(REC0826ch2(:,1),REC0826ch2(:,2), 'r');
% hold on
% plot(REC0826ch2(:,1),filtered_accel,'k');


% figure; plot(REC0826ch2(:,1),REC0826ch2(:,2), 'r');
% hold on
% plot(REC0826ch2(:,1),filtered_accel,'k');

figure; plot(REC0826ch2.REC0826ch2(:,1),REC0826ch2.REC0826ch2(:,2), 'r','DisplayName','raw acceleration');
title('Comparision between filted acceleration and raw accleration');
xlabel('Time (s)');
ylabel('Acceleration (m/s^2)');
hold on
plot(REC0826ch2.REC0826ch2(:,1),filtered_accel,'k','DisplayName','filtered acceleration');
hold off
legend('show');

REC0826ch1 = load('C:\Users\mm395\Documents\MATLAB\REC0826ch1.mat'); % load acceleration data
figure; plot(REC0826ch1.REC0826ch3(:,1),REC0826ch1.REC0826ch3(:,2), 'r','DisplayName','CoCo 80 Displacement');
hold on
plot(REC0826ch1.REC0826ch3(:,1),AD_V_acce,'k','DisplayName','Matlab Displacement');
xlabel('Time (s)');
ylabel('Displacement (m)');
hold off
legend('show');