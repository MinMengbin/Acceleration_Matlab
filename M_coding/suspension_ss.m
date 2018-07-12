%% this coding is used for simulation of a passive suspension system
% reference: http://ctms.engin.umich.edu/CTMS/index.php?example=Suspension&section=SystemModeling
clear all;

%% Set up
parameters_modelv3; % load parameters
% simutime = 92; % simulation time length
% simufre = 1/3200; %  simulation frequency
%% transfer function
s = tf('s');


%% G1 - G4 are the models for a secondary passive suspension system and
% spring tyre
% force control input (x1-x2) /U
G1 = ((M1+M2)*s^2+b2*s+K2)/((M1*s^2+b1*s+K1)*(M2*s^2+(b1+b2)*s+(K1+K2))-(b1*s+K1)*(b1*s+K1));

% force control input x1/U
G10 =  (M2*s^2+b2*s+K2)/((M1*s^2+b1*s+K1)*(M2*s^2+(b1+b2)*s+(K1+K2))-(b1*s+K1)*(b1*s+K1));

% % road disturbance input, output is the x1-x2  (x1-x2)/W
G2 = (-M1*b2*s^3-M1*K2*s^2)/((M1*s^2+b1*s+K1)*(M2*s^2+(b1+b2)*s+(K1+K2))-(b1*s+K1)*(b1*s+K1));
% 
% % road disturbance input, output is the sprung mass displacement x1/W
G3 =  (b1*b2*s^2+(b1*K2+b2*K1)*s+K1*K2)/((M1*s^2+b1*s+K1)*(M2*s^2+(b1+b2)*s+(K1+K2))-(b1*s+K1)*(b1*s+K1));
% % road disturbance input, output is the unsprung mass displacement x2/W
G4 =  (M1*b2*s^3+(M1*K2+b1*b2)*s^2+(b1*K2+b2*K1)*s+K1*K2)/((M1*s^2+b1*s+K1)*(M2*s^2+(b1+b2)*s+(K1+K2))-(b1*s+K1)*(b1*s+K1));



% %% Tyre model x/W  Frequency Response, Time response
%%% tyre model parameters
% kt = 100000; % tyre stiffness
% m = 5; % the mass of the hub
% G = kt/(m*s^2+2*kt); %% tyre model; 
% figure(1);
%  w = logspace(-1,5,10000);
% % frequency response of No suspension system
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
%     text(fs(1,1),20,['natural frequency is:' num2str(fs(1,1)/(2*pi),'%2.3f') 'Hz'] ,'Color','black','FontSize',14);   
%     subplot(3,1,2);
%     semilogx(w, squeeze(phase));
%     xlabel('Normalized Frequency rad/s'); % x-axis label
%     ylabel('Phase degree'); % y-axis label
%     subplot(3,1,3);
%     step(G);
% %% Tyre model x/W  Frequency Response, Time response




% %% No suspension system X1/W  Frequency Response, Time response
%  G5 = 1/(1+(M1+M2)*s^2/K2); 
% figure(1);
%  w = logspace(-1,3,10000);
% frequency response of No suspension system
% for i = 100:10:10000
%     G5 = 1/(1+(M1+M2)*s^2/i);   
%     [mag,phase] = bode(G5,w);
%     subplot(3,1,1);
%     plot(100:10:i,100:10:i);
%     title(['Current spring stiffness value is'  num2str(i,'%2.3f') 'N/m']);
%     xlabel('Number'); % x-axis label
%     ylabel('Value of Spring stiffness'); % y-axis label
%     subplot(3,1,2);
%     loglog(w, squeeze(mag));
%     title('Bode Diagram of No suspension system');
%     xlabel('Normalized Frequency rad/s'); % x-axis label
%     ylabel('Magnitutide dB'); % y-axis label
%     text(sqrt(i/(M1+M2)),20,['natural frequency is:' num2str(sqrt(i/(M1+M2)),'%2.3f') 'rad/s'] ,'Color','black','FontSize',14);   
%     subplot(3,1,3);
%     semilogx(w, squeeze(phase));
%     xlabel('Normalized Frequency rad/s'); % x-axis label
%     ylabel('Phase degree'); % y-axis label
%     drawnow
% end
% %% No suspension system X1/W  Frequency Response, Time response

% %% Passive suspension system X1/W  Frequency Response, Time Domain Analysis
    % % road disturbance input, x1/W
    G6 =  (b1*b2*s^2+(b1*K2+b2*K1)*s+K1*K2)/((M1*s^2+b1*s+K1)*(M2*s^2+(b1+b2)*s+(K1+K2))-(b1*s+K1)*(b1*s+K1));
    figure(3);
    w = logspace(-1,5,10000);
    % frequency response of passive suspension system
    for i = 100:100:10000
        b1 = i;
        G3 = (b1*b2*s^2+(b1*K2+b2*K1)*s+K1*K2)/((M1*s^2+b1*s+K1)*(M2*s^2+(b1+b2)*s+(K1+K2))-(b1*s+K1)*(b1*s+K1));  
        [mag,phase] = bode(G3,w);
        subplot(5,1,1);
        plot(100:10:i,100:10:i);
        % changing spring stiffness
    %     title(['Current spring stiffness value is '  num2str(i,'%2.3f') 'N/m']);
        % changing damper co-efficient
        title(['Current co-efficient of damper value is '  num2str(i,'%2.3f') 'N*s/m']);
        xlabel('Number'); % x-axis label
        ylabel('Value of Spring stiffness'); % y-axis label
        subplot(5,1,2);
        loglog(w, squeeze(mag));
        title('Bode Diagram of Passive suspension system');
        xlabel('Normalized Frequency rad/s'); % x-axis label
        ylabel('Magnitutide dB'); % y-axis label
       % text(sqrt(i/(M1+M2)),20,['Frequency is:  ' num2str(sqrt(i/(M1+M2)),'%2.3f') 'rad/s'] ,'Color','black','FontSize',14);   
        subplot(5,1,3);
        semilogx(w, squeeze(phase));
        xlabel('Normalized Frequency rad/s'); % x-axis label
        ylabel('Phase degree'); % y-axis label
        subplot(5,1,4);
        step(0.02*G3);
        title('Step Response of Passive Suspension System');
        xlabel('Time s'); % x-axis label
        ylabel('Displacement mm'); % y-axis label
        %subplot(5,1,5);
        drawnow
    end
    
%     % road disturbance input, x1/W Step response one figure
%     G6 =  (b1*b2*s^2+(b1*K2+b2*K1)*s+K1*K2)/((M1*s^2+b1*s+K1)*(M2*s^2+(b1+b2)*s+(K1+K2))-(b1*s+K1)*(b1*s+K1));
%     figure(3);
%     title('Step Response of Passive Suspension System');
%     xlabel('Time s'); % x-axis label
%     ylabel('Displacement mm'); % y-axis label
%     % frequency response of passive suspension system
%     for i = 1.0954e+03:100:1.0954e+03
%         b1 = i;
%         G3 = (b1*b2*s^2+(b1*K2+b2*K1)*s+K1*K2)/((M1*s^2+b1*s+K1)*(M2*s^2+(b1+b2)*s+(K1+K2))-(b1*s+K1)*(b1*s+K1));  
%         step(0.02*G3);
%         hold on
%         drawnow
%     end
    
    
    %%% Overshot and Setting Time at Different Damperness
%     figure(2);
%     comOT = animatedline;
%     view(3);
%     title('Overshot and Setting Time at Different Damperness');
%     xlabel('Value of Dampness'); % x-axis label
%     ylabel('Overshoot Value'); % y-axis label
%     zlabel('Settling Time s'); % y-axis label
%     for i = 100:100:5000
%         b1 = i;
%         G3 = (b1*b2*s^2+(b1*K2+b2*K1)*s+K1*K2)/((M1*s^2+b1*s+K1)*(M2*s^2+(b1+b2)*s+(K1+K2))-(b1*s+K1)*(b1*s+K1));  
%         stepinf = stepinfo(G3);
%         addpoints(comOT,i,stepinf.Overshoot,stepinf.SettlingTime);
%         drawnow limitrate
%     end
    %%% Overshot and Setting Time at Different Damperness


% %% Passive suspension system X1/W  Frequency Response, Time Domain Analysis


%% Skyhook suspension system X1/W  Frequency Response, Time Domain Analysis
% 
% G7 = k1*k2/((m1*s^2+b1*s+k1)*(m2*s^2+k1+k2)-k1*k1);
%     figure(4);
% %% step response at different dampness
%     title('Step Response of Skyhook Suspension System');
%     xlabel('Time s'); % x-axis label
%     ylabel('Displacement mm'); % y-axis label
% %% step response at different dampness
% 
% %% frequency response of Skyhook suspension system
%       w = logspace(-1,5,10000);
% %% frequency response of Skyhook suspension system
% 
%     for i = 100:10:1000
%         b1 = i;
%         G7 = k1*k2/((m1*s^2+b1*s+k1)*(m2*s^2+k1+k2)-k1*k1);
%         step(0.02*G7);  % step response Time Domain 
%         bode(G7,w);   % frequency response
%         legend(['Csky = ' num2str(i,'%2.f')]); 
%         hold on
%         drawnow
%     end
%     
    

% % %% Groundhook suspension system X1/W  Frequency Response, Time Domain Analysis
% G8 = k1*k2/((m2*s^2+b1*s+k1+k2)*(m1*s^2+k1)-k1*k1);
%    figure(5);
% %%% step response at different dampness
% %     title('Step Response of Groundhook Suspension System');
% %     xlabel('Time s'); % x-axis label
% %     ylabel('Displacement mm'); % y-axis label
% %%% step response at different dampness
% 
% %%% frequency response of Skyhook suspension system
%       w = logspace(-1,2,10000);
% %%% frequency response of Skyhook suspension system
% 
%     for i = 1000:400:1000
%         b1 = i;
%         G8 = k1*k2/((m2*s^2+b1*s+k1+k2)*(m1*s^2+k1)-k1*k1);
% %         step(0.02*G8);  % step response Time Domain 
%         bode(G8,w);   % frequency response
% %         legend(['Csky = ' num2str(i,'%2.f')]); 
%         hold on
%         drawnow
%     end
    
    
% X1/W Comparison Between No Suspension, Passive, Skyhook, Groundhook 
% % both Frequency Response and Time Domain)
G5 = 1/(1+M1*s^2/K2);  % Model for No suspension
% G3 Model for Passive suspension
G7 = k1*k2/((m1*s^2+b1*s+k1)*(m2*s^2+k1+k2)-k1*k1); % Model for Skyhook 
G8 = k1*k2/((m2*s^2+b1*s+k1+k2)*(m1*s^2+k1)-k1*k1); % Model for Groundhook



% X2/W Comparison Between No Suspension, Passive, Skyhook, Groundhook 
% % both Frequency Response and Time Domain)

%%% Model for Passive suspension
% G4 Model for Passive suspension
G6 = (m1*s^2+b1*s+k1)*k2/((m1*s^2+b1*s+k1)*(m2*s^2+k1+k2)-k1*k1); % Model for Skyhook 
G9 = (m1*s^2+k1)*k2/((m2*s^2+b1*s+k1+k2)*(m1*s^2+k1)-k1*k1); % Model for Groundhook



figure;
% w = logspace(-1,2,100);
w = (0.01:0.001:50*2*pi);
 
% options = bodeoptions;
% options.FreqUnits = 'Hz'; % or 'rad/second', 'rpm', etc.
% frequency response of No suspension system

%%bode diagram
% bode(G5,w);
% hold on
% bode(G3,w);
% hold on
% bode(G7,w);
% hold on
% bode(G8,w);
% hold on
% legend('No Suspesion','Passive Suspension','Ideal Skyhook', 'Ideal Groundhook');
% hold off
%bode diagram
fs = (0.5/pi)*sqrt(kt/m);

%% step reponse comparison
figure;
% step(G5);
% hold on
% step(G3);
% hold on
% step(G7);
% hold on
% step(G8);
% hold on
% legend('No Suspesion','Passive Suspension','Ideal Skyhook', 'Ideal Groundhook');
% hold off


% frequency analysis for x2/w
% figure;
% bode(G4);
% hold on
% bode(G6);
% hold on
% bode(G9);
% hold on
% legend('Passive Suspension','Ideal Skyhook', 'Ideal Groundhook');
% hold off

w1 = sqrt(k1/m1);
f1 = (0.5/pi)*sqrt((k1+k2)/m1);
w2 = sqrt(k2/m2);
f2 = (0.5/pi)*sqrt(k2/m2);
% ModelS = [G5 G6 G7 G8];
% 
% for i = 1:1:4
% [mag,phase] = bode(ModelS(i),w);
% subplot(2,1,1);
% plot(w/(2*pi), squeeze(mag));
% hold on
% title('Bode Diagram of Passive suspension system');
% xlabel('Normalized Frequency Hz'); % x-axis label
% ylabel('|X1/Y|'); % y-axis label
% %text(sqrt(i/(M1+M2)),20,['Frequency is:  ' num2str(sqrt(i/(M1+M2)),'%2.3f') 'rad/s'] ,'Color','black','FontSize',14);   
% subplot(2,1,2);
% plot(w/(2*pi), squeeze(phase));
% xlabel('Normalized Frequency Hz'); % x-axis label
% ylabel('Phase degree'); % y-axis label
% legend('No Suspesion','Passive Suspension','Skyhook', 'Groundhook');
% hold on
% end
% hold off
%% X1/W Comparison Between No Suspension, Passive, Skyhook, Groundhook 


%% open loop step response analysis
% step(G1); % force input
% t = (0:0.01:simutime);
% stepd = frest.createStep('Ts',simufre,'StepTime',0.1,'StepSize',1,'FinalTime',simutime);
% u = stepd.data; % 
%  y = step(G2,t); % gound disturbance input with magnitude 0.1m
% u = rand(1, simutime/simufre+1); % random input

% gound disturbance input of field collected data
%  groundinput = load('groundinput.mat'); % load acceleration data
%  u = 100 * groundinput.ans(2,:);
%  t = groundinput.ans(1,:);
%  yx1 = lsim(G3,u,t); % road disturbance response
%  yx2 = lsim(G4,u,t); % road disturbance response
%  y = step(G2, t); % step response 
%  figure; 
%  plot(t,y,':.b');
%  hold on;
%  passivesus = load('passivesus.mat'); % load acceleration data
%  plot(passivesus.ans(1,:),passivesus.ans(2,:),'--.g');
%  hold off;


% animate response

% z0 = u;                     % road elevation
% z1 = yx2;                   % unspring mass cm position
% z2 = yx1;                   % sprung mass position
% zmf = 1;    % exaggerate response for better visualization
% for i=1:500:length(t)
%     passivets([z0(i), z1(i)*zmf, z2(i)*zmf, t(i)],t,u,t(i),1,simutime);
%     text(-simutime/2,7,['Ponit number:' num2str(i,'%2.f')],'Color','black','FontSize',14);
%     drawnow
% end

% figure(2);
% plot(t,yx1,':.r');
% hold on;
% plot(t,u,'--.g');
% title('Compare Between Road Disturbance and Suspended Mass Displacement')
% xlabel('Time /s'); % x-axis label
% ylabel('Displacement /cm'); % y-axis label
% legend('Suspended Mass','Road Disturbance')
% hold off;

