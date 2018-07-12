% rxCh = canChannel('Kvaser', 'Leaf Light v2 1', 1);
% configBusSpeed(rxCh,250000);
% get(rxCh);
start(rxCh);
% rxMsg = receive(rxCh, Inf, 100);
% rxMsg(1:25, :);
% Create figure
pause (0.1);
Accel_z = figure;

% % Create axes
% axes1 = axes('Parent',Accel_z);
% hold(axes1,'on');

hold on;
% Create plot
%plot(X1,Y1,'Marker','diamond','LineWidth',2);

% Create xlabel
xlabel('Timestamp (s)');

% Create ylabel
ylabel('Accel (m/s/s)');

% Create title
title('Accel-z vs Time');

%box(axes1,'on');
box on;
plot(0,0,'-or')
j = 0;
for c = 1: 100
    message = receive(rxCh,1);
    while (size(message) == 0)
       message = receive(rxCh,1);
    end
 %   if message(1,1).ID == 486474880
%         if j == 0
%             firstT = message(1,1).Timestamp;  % first one time stamp
%             j = j+1;
%         end
       %plot(message(1,1).Timestamp - firstT,message(1,1).Timestamp - firstT,'-or')
        plot(message(1,1).Timestamp,message(1,1).Timestamp,'-or')
        %fprintf('timestamp is %d\n',message(1,1).Timestamp);
        pause(0.1); %time delay for showing 
%    end
%     for i = 1: 10 
%     end
%     plot(message(1,1).Timestamp,message(1,1).Data(1,1),'-or')
end
% for c = 1 : 20
   
% end
stop(rxCh);