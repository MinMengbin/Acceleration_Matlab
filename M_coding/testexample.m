% %%%function example x = x+6. 
% % function x = testexample(x,~)
% % x = x+6;
% % %disp(x);
% % end
% %%%function example x = x+6.  
% 
% %%%ODE example dx/dt = y.  
% delta = 1;
% F = inline('y','t','y');
% %formula(F);
% opts = odeset('RelTol',1.e-4);
% ode45(F,[0 2/delta],delta,opts);
% %%%ODE example dx/dt = y. 
% 
% %%%ODE example dx/dt = -3 exp(-t). 
% % function testexample
% % 
% % % SOLVE  dx/dt = -3 exp(-t).  
% % % initial conditions: x(0) = 0
% % 
% % t=0:0.01:10;   % time scalex
% % initial_x=0;
% % 
% % [t,x]=ode45( @rhs, t, initial_x);
% % 
% % plot(t,x);
% % xlabel('t'); ylabel('x');
% % 
% % function dxdt=rhs(t,x)
% %         dxdt = 3*exp(-t);
% %     end
% % end
% %%%%ODE example dx/dt = -3 exp(-t). 




%% test linmode
clear all;

%% Set up
% parameters_modelv3; % load parameters
PIDcontrol;
[A,B,C,D] = linmod('fullyactive');
% [num,den] = ss2tf(A,B,C,D);
sys = ss(A,B,C,D)
 size(sys)
%[num1,den1] = ssdata(G6);