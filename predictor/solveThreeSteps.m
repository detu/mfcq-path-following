function [x_init, y_init, qp_exit, elapsedqp, deltaT] = solveThreeSteps(prob, p, x_init, y_init, step, lb, ub, N, x0, delta_t, oldEta)
%SOLVETHREESTEPS Summary of this function goes here
% 
% Solve 3 step described in MFCQ paper  
%
% [OUTPUTARGS] = SOLVETHREESTEPS(INPUTARGS) Explain usage here
% 
% Examples: 
% 
% Provide sample usage code here
% 
% See also: List related files here

% $Author: suwartad $	$Date: 2017/05/11 14:24:20 $	$Revision: 0.1 $
% Copyright: Process Control Group - NTNU Trondheim 2017

% TO DO:
% - compute the Hessian of constraints (Hc) ONE TIME ONLY 
% - TRY compute Hc only at the BEGINNING OF MPC iteration

% Derivatis w.r.t to deltaT can be computed ONE TIME.
global flagDt mpciter;

if mpciter == 1
    flagDt = 1;
else
    flagDt = 0;
end
% obtain derivatives information
[~,g,H,Lxp,cin,~,~,Jeq,dpe,~,Hc] = prob.obj(x_init,y_init,p, N);


% First: Corrector Step
[deltaXc,deltaYplus]      = solveCorrectStep(H, Jeq, g, cin, y_init);

% Second: Predictor
[deltaXp,deltaYp, exitQP] = solveQPPredict(H, Jeq, g, cin, Hc, step, dpe, lb, ub, deltaXc);

% Sum up steps from Corrector and Predictor
deltaX = deltaXc + deltaXp;
deltaY = deltaYplus + deltaYp;

% Calculate new Eta 
xCurrent    = x_init + deltaX;
yCurrent    = y_init.lam_g + deltaY;
flagDt      = 0;
[~,g,~,~,cin,~,~,Jeq,~,~,~] = prob.obj(xCurrent,yCurrent,p, N);
[newEta, z] = computeEta(Jeq, g, yCurrent, cin);

% Update deltaT  
deltaT      = updateDeltaT(oldEta, newEta, delta_t);

% Need not to update active-set since we have equality constraint

% Third: Jump Step
[y_init,exitLP]   = solveJumpLP(Jeq, Lxp, g, dpe, cin, y_init, step, z);


end

function [dXc,dYplus] = solveCorrectStep(H, Jeq, g, c, y)

[n, m] = size(Jeq);
lhs     = [H                   -Jeq'; ...
           Jeq                 zeros(n,n)];
rhs     = -[g-Jeq'*y.lam_g ; c];
solCorrectStep = lhs\rhs;
dXc    = solCorrectStep(1:m);
dYplus = solCorrectStep(m+1:end);

end

function [dXp, dYp, qp_exit] = solveQPPredict(H, Jeq, g, cin, Hc, step, dpe, lb, ub, deltaXc)
% QP consists of equality and bound constraints 

% QP setup
% A   = [];
% b   = [];
f   = [];

% Only Equality Constraint
ceq         = cin;
[numY,numX] = size(Jeq);
Hcx         = zeros(numY,numX);
for i=1:numY
    Hcx(i,:) = (Hc{:,:,i}*deltaXc)';
end
Aeq  = Jeq + Hcx;
beq  = dpe*step + ceq;   %OK

% TOMLAB setup
Prob   = qpAssign(H, f, Aeq, beq, beq, lb, ub);
Prob.optParam.eps_x = 1e-7;
Prob.optParam.fTol  = 1e-7;
Prob.optParam.xTol  = 1e-7;
startqp  = tic;
Result = tomRun('qp-minos', Prob, 1);

elapsedqp = toc(startqp);
fprintf('QP solver runtime: %f\n',elapsedqp);
qp_exit = Result.ExitFlag;
if qp_exit == 0
    dXp       = Result.x_k;
    %qp_val  = Result.f_k;
else
    keyboard;
end

% numX        = size(x_init,1);
% lamda.lam_x = Result.v_k(1:numX);
% lamda.lam_g = -Result.v_k(numX+1:end);
% fprintf('QP return: %d\n', qp_exit);
dYp = -Result.v_k(numX+1:end);

% % CPLEX setup
% options = cplexoptimset;
% options.Display = 'on';
% startqp  = tic;
% [dXp, fval, exitflag, output, lambda] = cplexqp (H, [], [ ], [ ], Aeq, beq, lb, ub, [ ], options);
% 
% elapsedqp = toc(startqp);
% fprintf('QP solver runtime: %f\n',elapsedqp);
% fprintf ('\nSolution status = %s \n', output.cplexstatusstring);
% fprintf ('Solution value = %f \n', fval);

end

function [lpSol, exitflag] = solveJumpLP(Jeq, Lxp, g, dpe, cin, y, step, z)
% test LP with GUROBI and CLPEX


numY = size(Jeq,1);

lb       = -Inf*ones(numY,1);
ub       = Inf*ones(numY,1);
f        = dpe*step;
A        = [Jeq';-Jeq'];
b        = [g+abs(z);abs(z)-g];
%option   = optimoptions('linprog','Algorithm','dual-simplex','Display','off');
option = cplexoptimset;
option.Display = 'on';
startlp  = tic;
%[lpSol,~,exitflag] = linprog(f, A, b, [], [], lb, ub, y.lam_g, option );
%[lpSol,~,exitflag] = linprog(f, A, b, [], [], lb, ub, [], option );

[lpSol,~,exitflag] = cplexlp(f, A, b, [], [], lb, ub, [], option );

elapsedlp = toc(startlp);
fprintf('LP solver runtime: %f\n',elapsedlp);

end

function deltaT = updateDeltaT(currentEta, nextEta, t)
gamma = 0.7;
alpha = 0.6;

% check good eta condition, if satisfied increase delta_t
if (nextEta < currentEta^(1+gamma))
    %if (nextEta < currentEta)
    deltaT = min(1-t,t/alpha);
else
    deltaT = min(1-t,t);
end

end
