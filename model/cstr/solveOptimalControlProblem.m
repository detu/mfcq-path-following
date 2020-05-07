function [u, lamda, lbw, ubw, objVal, params, elapsednlp] = solveOptimalControlProblem(optProblem, N, x0, u0, x0_measure)
%SOLVEOPTIMALCONTROLPROBLEM Summary of this function goes here
% 
% [OUTPUTARGS] = SOLVEOPTIMALCONTROLPROBLEM(INPUTARGS) Explain usage here
% 
% Examples: 
% 
% Provide sample usage code here
% 
% See also: List related files here

% $Author: suwartad $	$Date: 2018/04/19 04:28:41 $	$Revision: 0.1 $
% Copyright: Process Control Group - NTNU Trondheim 2018

import casadi.*
% call ODE15s N-times initial guess in optimization
x(1,:) = x0';
for k=1:N
    x(k+1,:) = x0';    % ONE-TIME SIMULATION !
end

[J,g,w0,w,lbg,ubg,lbw,ubw,params] = optProblem(x, u0, N, x0_measure);
prob = struct('f', J, 'x', vertcat(w{:}), 'g', vertcat(g{:}));
options = struct;
options.ipopt.tol                = 1e-12;
options.ipopt.constr_viol_tol    = 1e-12;
solver = nlpsol('solver', 'ipopt', prob, options);

% Solve the NLP
startnlp   = tic;
sol        = solver('x0', w0, 'lbx', lbw, 'ubx', ubw, 'lbg', lbg, 'ubg', ubg);
elapsednlp = toc(startnlp);
fprintf('IPOPT solver runtime = %f\n',elapsednlp);
success = strcmp(solver.stats.return_status,'Infeasible_Problem_Detected');
if (success)
    keyboard;
end

u      = full(sol.x);
lamda  = full(sol.lam_g);
objVal = full(sol.f);

end
