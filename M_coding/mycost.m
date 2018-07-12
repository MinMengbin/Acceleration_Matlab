function cost = mycost(k)
% This function takes two parameters k = [k1 k2] in the system and computes
% the cost funtion.

% Simulation time is defined as 10s
tf = 10;

% Initial condition is ARBITRRY defined as x0 = [0.1 -0.2];
x0 = [0.1 -0.2];

% Simulate the system with different k1 and k2
[~, y] = ode45(@(t,x) mysys(t,x,k(1),k(2)), [0 tf], [x0 0]);

% The cost function is the WEIGHTED SUM OF THE INTEGRATED NORM OF x(1) and x(2)
cost = y(end);
end