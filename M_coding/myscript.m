% Initial guess for the parameters
k1 = 1;
k2 = 1;

% Specify the lower and upper bounds
klow = [0.1 0.1];
kup  = [10   10];

% Solve the optimizaiton proble mx = fmincon(fun,x0,A,b,Aeq,beq,lb,ub)
opt = optimoptions('fmincon');
opt.Display = 'Iter';
optk = fmincon(@mycost, [k1 k2], [], [], [], [], klow, kup, [],opt);

% Display results
disp('The optimal values are')
disp(optk)