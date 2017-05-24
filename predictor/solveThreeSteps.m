function [x_init, y_init, qp_exit, deltaT, newEta, success, p_t] = solveThreeSteps(prob, p, p_final, x_init, y_init, step, lb, ub, lb_init, ub_init, N, x0, t, delta_t, p_0, p_t, oldEta)
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
% - READ active bound constraint as equality constraint !

% Derivatis w.r.t to deltaT can be computed ONE TIME.
global flagDt mpciter;

if ((mpciter == 1) && (t == 0))
    flagDt = 1;
else
    flagDt = 0;
    %load Hc.mat;
end
% obtain derivatives information
[~,g,H,Lxp,cin,~,~,Jeq,dpe,~,Hc] = prob.obj(x_init,y_init,p, N);
if isempty(Hc)
    load Hc.mat;
end


% First: Corrector Step
[deltaXc,deltaYplus]      = solveCorrectStep(H, Jeq, g, cin, y_init);

% Second: Predictor
[deltaXp,deltaYp, exitQP] = solveQPPredict(H, Jeq, g, cin, Hc, step, dpe, lb, ub, deltaXc);

% Sum up steps from Corrector and Predictor
deltaX = deltaXc + deltaXp;
deltaY = deltaYplus + deltaYp.lam_g;

% Calculate new Eta 
xCurrent       = x_init + deltaX;
yCurrent.lam_g = y_init.lam_g + deltaY;
yCurrent.lam_x = y_init.lam_x + deltaYp.lam_x;
%yCurrent.lam_x = y_init.lam_x;
% yCurrent.lam_g = deltaY;
% yCurrent.lam_x = deltaYp.lam_x;
flagDt         = 0;
[~,g,~,~,cin,~,~,Jeq,~,~,~] = prob.obj(xCurrent,yCurrent,p, N);
[newEta, z] = computeEta(Jeq, g, yCurrent, cin);

% Checking condition (5.1) in the paper !!!
success = 0;
while (newEta > max(0.2,oldEta))  % EtaMax = 1e-2, should be 1e-6?
    
    % decrease step
    delta_t = 0.6*delta_t;
    tk      = t + delta_t;
    p_t     = (1 - tk)*p_0 + tk*p_final;
    step    = p_t - p;
    
    % update bound constraint
    if(~isempty(lb_init))
        lb = lb_init - deltaX; 
        ub = ub_init - deltaX; 
    else
        lb = [];
        ub = [];
    end
    
    % solve Predictor
    [deltaXp,deltaYp, exitQP] = solveQPPredict(H, Jeq, g, cin, Hc, step, dpe, lb, ub, deltaXc);
    
    % Sum up steps from Corrector and Predictor
    deltaX = deltaXc + deltaXp;
    deltaY = deltaYplus + deltaYp.lam_g;
    
    % Calculate new Eta
    xCurrent       = x_init + deltaX;
    yCurrent.lam_g = y_init.lam_g + deltaY;
    yCurrent.lam_x = y_init.lam_x + deltaYp.lam_x;
    %yCurrent.lam_g = deltaY;
    flagDt         = 0;
    [~,g,~,~,cin,~,~,Jeq,~,~,~] = prob.obj(xCurrent,yCurrent,p, N);
    [newEta, z] = computeEta(Jeq, g, yCurrent, cin);
end
success = 1;

% Update deltaT  
deltaT      = updateDeltaT(oldEta, newEta, delta_t);

% Need not to update active-set since we have equality constraint

% Third: Jump Step
%[y_init.lam_g,exitLP]   = solveJumpLP(Jeq, Lxp, g, dpe, cin, y_init, step, z);
[y_init.lam_g,exitLP]   = solveJumpLP(Jeq, Lxp, g, dpe, cin, yCurrent, step, z);

% update y_init.lam_x
y_init.lam_x = yCurrent.lam_x; % Is it Correct? 
% update y_init.lam_x just for the active bound constraint... 

% set dummy variables for time being
qp_exit   = 1;
%elapsedqp = 1;
x_init    = xCurrent;

end

function [dXc,dYplus] = solveCorrectStep(H, Jeq, g, c, y)

[n, m] = size(Jeq);
lhs     = [H                   -Jeq'; ...
           Jeq                 zeros(n,n)];
%rhs     = -[g-Jeq'*y.lam_g ; c];
rhs     = -[g + Jeq'*y.lam_g + y.lam_x; c];
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
%Aeq  = Jeq;
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
%dYp = -Result.v_k(numX+1:end);
dYp.lam_x = -Result.v_k(1:numX);
dYp.lam_g = -Result.v_k(numX+1:end);


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
% INCLUDE ACTIVE BOUND CONSTRAINT HERE !

%numY = size(Jeq,1);
[numY, numX] = size(Jeq);

% % active bound constraints
% boundMultiplier    = y.lam_x;
% positiveBoundMult  = abs(boundMultiplier);
% activeIndex        = find(positiveBoundMult>1e-3);
% numActiveBoundCons = numel(activeIndex);
% % Jbc                = eye(numActiveBoundCons);   
% % JbcExt             = zeros(numActiveBoundCons, numX-numActiveBoundCons);
% % JAbc               = [Jbc JbcExt];
% % JActive            = [Jeq;JAbc];

% JAbc = zeros(numActiveBoundCons,numX);
% for i=1:numActiveBoundCons
%     row         = activeIndex(i);
%     JAbc(i,row) = 1;
% end
% JActive         = [Jeq;JAbc];

lb       = -Inf*ones(numY,1);
ub       = Inf*ones(numY,1);
% lb       = -Inf*ones(numY+numActiveBoundCons,1);
% ub       = Inf*ones(numY+numActiveBoundCons,1);
f        = dpe*step;
% f        = [f;zeros(numActiveBoundCons,1)];
A        = [Jeq';-Jeq'];
% A        = [JActive';-JActive'];
%b        = [g+abs(z);abs(z)-g];
%b        = [abs(z)-g;abs(z)+g];
b        = [abs(z) - g - y.lam_x; abs(z) + g + y.lam_x];
% A        = [Jeq';-Jeq'];
%b        = [abs(z) + g + y.lam_x; abs(z) - g - y.lam_x];
% option   = optimoptions('linprog','Algorithm','dual-simplex','Display','off');
option = cplexoptimset;
option.Display = 'on';
startlp  = tic;
%[lpSol,~,exitflag] = linprog(f, A, b, [], [], lb, ub, y.lam_g, option );
%[lpSol,~,exitflag] = linprog(f, A, b, [], [], lb, ub, [], option );

[lpSol,~,exitflag] = cplexlp(f, A, b, [], [], lb, ub, [], option );

elapsedlp = toc(startlp);
% lpSol     = -lpSol; % different sign !
fprintf('LP solver runtime: %f\n',elapsedlp);

% Take out active bound constraint here !


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
