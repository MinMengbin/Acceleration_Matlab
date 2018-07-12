% Define the STATE-SPACE representation of a simple 2nd order system
% The system is 
% dx1/dt = x2;
% dx2/dt = -k1*x1 - k2 * x2;
function dx = mysys(~,x,k1,k2)
dx = zeros(3, 1); % Initialization of the derivative
dx(1) = x(2);
dx(2) = -k1 * x(1) - k2 * x(2);
% The augmented state x3 is for the purpose of doing the weighted
% sum of the norm of x1 and x2.
% This cost function is very similar, if not the same, as in LQR.
% The weights I chose are simpley 1 and 2.
dx(3) = 1*x(1)*x(1) + 2*x(2)*x(2);
end