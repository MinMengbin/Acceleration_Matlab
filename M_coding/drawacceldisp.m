% this piece of coding is to plot the average peak value of acceleration &
% displacement at different speeds
speedlist = [0.4, 0.6, 1.0, 1.7]; % 0 engin off; 0 engin on; 0.4 m/s
peakacceleration = [1.5, 3.4, 6.1, 11.2];% m/s^2
peakdisplacement = [0.004, 0.005, 0.006, 0.008]; % m

EngineS = ["OFF", "ON"];
accelengine = [0.02, 1];% m/s^2
dispengine = [0.000, 0.001]; % m

% plot the peak acceleration value relative to speed
figure; plot(speedlist, peakacceleration,'-.*B');
title('Peak value of acceleration at different speeds');
xlabel('Speed (m/s)');
ylabel('Acceleration (m/s^2)');
legend('Peak Acceleration Value','Location','northeast')

% plot the peak displacement value relative to speed
figure;plot(speedlist, peakdisplacement,'-.*r');
title('Peak value of diplacement at different speeds');
xlabel('Speed (m/s)');
ylabel('Displacement (m)');
legend('Peak Acceleration Displacement','Location','northeast')

% plot the peak acceleration value relative to engine status
figure; plot(accelengine(1),'-.*B');
set(gca,'xticklabel',EngineS(1));
title('Peak value of acceleration at different engine status');
xlabel('Engine Status');
ylabel('Acceleration (m/s^2)');
hold on
plot(accelengine(2),'-.*B');
set(gca,'xticklabel',EngineS(2));
hold off

% plot the peakacceleration value relative to engine status
figure; plot(dispengine,'-.*B');
set(gca,'xticklabel',EngineS);
title('Peak value of diplacement at different engine status');
xlabel('Engine Status');
ylabel('Acceleration (m/s^2)');