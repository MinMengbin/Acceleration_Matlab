% state-space representation of passive susension system
% dx(t)/dt = Ax(t) + Bu(t)
% y(t) = Cx(t) + Du(t)

%        |  xm(t)      |  sprung mass displacement
% x(t) = |  dxm(t)/dt  |  sprung mass velocity
%        |  xmu(t)     |  unsprung mass displacement
%        |  dxmu(t)/dt |  unsprung mass velocity

%initialise system parameters
k_t = 100000; %  the stiffness of the tyre is 100000 n/m
m = 150; % frame mass is 150 kg
m_u = 30; % unprung frame mass is 30 kg
%ks; %  the constant of spring
%dc; %  the coefficient of the damper

function dx = syspassive(~,x,x_in,k_s,d_c)
disp(' This is a state-space represenstation for a passive suspension system')

%dx = zeros(4,1);

A = [   0         1              0             0;
     -k_s/m     -d_c/m          k_s/m         d_c/m;
        0        0               0              1;
     k_s/m_u     d_c/m_u    -(k_s+k_t)/m_u  -d_c/m_u ];


B = [      0;
           0;
           0;
        k_t/m_u];
   
dx = dot(A,x)+ dot(B,x_in);



end