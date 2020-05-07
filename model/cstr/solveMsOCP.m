function [u, lamda, lbw, ubw, objVal, params] = solveMsOCP(optProblem, system, N, t0, x0, u0, T, mpciter, u_nlp, x_nlp, x0_measure)
%SOLVEMSOCP Summary of this function goes here
% 
% [OUTPUTARGS] = SOLVEMSOCP(INPUTARGS) Explain usage here
% 
% Examples: 
% 
% Provide sample usage code here
% 
% See also: List related files here

% $Author: suwartad $	$Date: 2018/05/04 05:31:49 $	$Revision: 0.1 $
% Copyright: Process Control Group - NTNU Trondheim 2018

import casadi.*

% call ODE15s N-times initial guess in optimization
x(1,:) = x0';
for k=1:N
    x(k+1,:) = x0';    % ONE-TIME SIMULATION !
end

[J,g,w0,w,lbg,ubg,lbw,ubw,params] = optProblem(x, u0, N, x0_measure);

prob    = struct('f', J, 'x', vertcat(w{:}), 'g', vertcat(g{:}));
options = struct;
options.ipopt.tol                = 1e-12;
options.ipopt.constr_viol_tol    = 1e-10;
solver = nlpsol('solver', 'ipopt', prob, options);

% Solve the NLP
startnlp = tic;
sol   = solver('x0', w0, 'lbx', lbw, 'ubx', ubw, 'lbg', lbg, 'ubg', ubg);
elapsednlp = toc(startnlp);
fprintf('IPOPT solver runtime = %f\n',elapsednlp);

u           = full(sol.x);
lamda.lam_g = full(sol.lam_g);
lamda.lam_x = full(sol.lam_x);
objVal      = full(sol.f);

end
