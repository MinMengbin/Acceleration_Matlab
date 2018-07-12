LF = 1.2;  % front hub displacement from body gravity center (m)
LR = 1.2;  % rear hub displacement from body gravity center (m)

b2 = 1000; % front suspension damping in (N sec/m)
b1 = 1000; % rear suspension damping in (N sec/m)

k2 = 16000; % front suspension stiffness in (N/m)
k1 = 16000; % rear suspension stiffness in (N/m)

% https://uk.mathworks.com/help/robust/gs/active-suspension-control-design.html
% Pnumatic tyre stiffness
k4 = 190000; % front tyre stiffness in (N/m) 
k3 = 190000; % rear tyre stiffness in (N/m)
cantoo
M_sb = 400;  % sprung body mass (kg)
M_UsbF = 30; % unsprung front body mass (kg)
M_UsbR = 30; % unsprung rear body mass (kg)
I = 188.76;    % body moment of inertia about y-axis in (kg m^2)
aspa_4w;