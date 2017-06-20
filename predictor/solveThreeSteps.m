function [x_init, y_init, qp_exit, delta_t, success] = solveThreeSteps(prob, x_init, y_init, step, lb, ub, N, x0, t, delta_t, p_0, p_t, oldEta, ub_init)
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
etamax  = oldEta;
success = 0;
%if ((mpciter == 1) && (t == 0))
if (t == 0)
    flagDt = 1;
else
    flagDt = 0;
    %load Hc.mat;
end
numX = size(x_init,1);
numY = size(y_init.lam_g,1);

%if oldEta > 1
if oldEta > 5
    % decrease deltaT
    delta_t = 0.6*delta_t;
    qp_exit = 0;
    % debug 
    fprintf("keyboard debug \n");
    keyboard;
else
    % proceed with everything... 
    
    % obtain derivatives information
    [~,g,H,Lxp,cin,~,~,Jeq,dpe,~,Hc] = prob.obj(x_init,y_init,p_0, N);
    if isempty(Hc)
        load Hc.mat;
    end
    
    % First: Corrector Step
    [deltaXc,deltaYplus,activeBoundInd,JAbc]      = solveCorrectStep(H, Jeq, g, cin, y_init, x_init, ub_init);
    
    % Second: Predictor
    [deltaXp,deltaYp, qp_exit] = solveQPPredict(H, Jeq, g, cin, Hc, step, dpe, lb, ub, deltaXc, activeBoundInd);
    
    % Calculate new Eta 
    xCurrent         = x_init + deltaXc + deltaXp;
    yCurrent.lam_g   = y_init.lam_g + deltaYplus.lam_g + deltaYp.lam_g;
    

    % update the dual's bound constraints
    %yCurrent.lam_x = y_init.lam_x + deltaYp.lam_x;
    yCurrent.lam_x = y_init.lam_x + deltaYp.lam_x + deltaYplus.lam_x;
    
    flagDt         = 0;
    [~,g,~,~,cin,~,~,Jeq,~,~,~] = prob.obj(xCurrent,yCurrent,p_t, N);
    [newEta, z] = computeEta(Jeq, g, yCurrent, cin);
    
    % Checking condition (5.1)
    if newEta <= max(5,etamax)   % Parameters.etaMax = 5.0
        % Update primal and dual variables 
        deltaX = deltaXc + deltaXp;
        deltaY = deltaYplus.lam_g + deltaYp.lam_g;
        
        % Update deltaT  
        delta_t = updateDeltaT(oldEta, newEta, delta_t);
        
        % Third: Jump Step
        [lpSol, exitLP]   = solveJumpLP(Jeq, Lxp, g, dpe, cin, yCurrent, step, z, JAbc);
        if exitLP >= 0
            y_init.lam_g = lpSol(1:numY);
            %y_init.lam_x = yCurrent.lam_x; % Is it Correct? 
            y_init.lam_x = zeros(numX,1);
            y_init.lam_x(activeBoundInd) = lpSol(numY+1:end);
        end
        
        success = 1;
        x_init  = x_init + deltaX;
        qp_exit = 1;
        
        % Update etaMax
        if newEta > etamax
            etamax = newEta;
        end
        
    else
        % decrease deltaT
        delta_t = 0.6*delta_t;
        qp_exit = 0;
        % debug
        fprintf("keyboard debug \n");
        keyboard;
    end
end

end

function [dXc,dYplus,activeIndex,JAbc] = solveCorrectStep(H, Jeq, g, c, y, x, ub)

% active bound constraints
boundMultiplier    = y.lam_x;
positiveBoundMult  = abs(boundMultiplier);
activeIndex        = find(positiveBoundMult>1e-1);  % set active bound constraint.
paramIndex         = find(activeIndex <= 84);       % remove the first 84 constraints (the parameter)
activeIndex(paramIndex) = [];
numActiveBoundCons = numel(activeIndex);

[numY,m] = size(Jeq);
JAbc     = zeros(numActiveBoundCons,m);
for i=1:numActiveBoundCons
    row         = activeIndex(i);
    JAbc(i,row) = 1;
end

Jeq     = [Jeq;JAbc];
n       = size(Jeq,1);
lhs     = [H                   -Jeq'; ...
           Jeq                 zeros(n,n)];
%rhs     = -[g-Jeq'*y.lam_g ; c];
%rhs     = -[g + Jeq'*y.lam_g + y.lam_x; c];
xActive = x(activeIndex);
rhs     = -[g + Jeq'*([y.lam_g; y.lam_x(activeIndex)]); [c;(ub(activeIndex)-xActive)]];
solCorrectStep = lhs\rhs;
dXc    = solCorrectStep(1:m);
dYplus.lam_g = solCorrectStep(m+1:m+numY);
dYplus.lam_x = zeros(m,1);
dYplus.lam_x(activeIndex) = solCorrectStep(m+numY+1:end);

end

function [dXp, dYp, elapsedqp] = solveQPPredict(H, Jeq, g, cin, Hc, step, dpe, lb, ub, deltaXc, activeBoundInd)
% QP consists of equality and bound constraints 

% TO DO: include active-bound constrain here !

% QP setup
A   = [];
b   = [];
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

% build equality constraint from active bound constraints
numBaC = size(activeBoundInd,1);
for i=1:numBaC
    % put strongly active constraint on boundary
    indB         = activeBoundInd(i);
    ub(indB)     = 0;        % keep upper bound on boundary
    lb(indB)     = 0;
end

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

% % QUADPROG setup
% %option = optimset('Display','off','Algorithm','active-set');
% option = optimset('Display','iter','Algorithm','interior-point-convex');
% startqp  = tic;
% [dXp, ~,  exitflag, ~, lambda] = quadprog(H, f, A, b, Aeq, -beq, lb, ub, [], option);
% elapsedqp = toc(startqp);
% fprintf('QP solver runtime: %f\n',elapsedqp);
% dYp.lam_g = -lambda.eqlin;
% dYp.lam_x = -lambda.lower - lambda.upper;

% % QPC setup
% dsp = 1;
% startqp  = tic;
% [dXp, error, lambda] = qpip(full(H), zeros(numX,1), [], [], Aeq, -beq, lb, ub, dsp);
% elapsedqp = toc(startqp);
% fprintf('QP solver runtime: %f\n',elapsedqp);
% dYp.lam_g = -lambda.eqlin;
% dYp.lam_x = -lambda.lower - lambda.upper;


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

% % QPOASES setup
% import casadi.*
% x      = MX.sym('x',numX);
% f_gn   = 0.5*(x'*(H)'*H*x);
% g_l    = beq + Aeq*x;
% qp     = struct('x',x, 'f',f_gn,'g',g_l);
% solver = qpsol('solver', 'qpoases', qp);
% startqp   = tic;
% sol       = solver('lbg',zeros(numY,1),'ubg',zeros(numY,1),'lbx',lb,'ubx',ub);
% elapsedqp = toc(startqp);
% fprintf('QP solver runtime: %f\n',elapsedqp);
% dXp       = full(sol.x);
end

function [lpSol, exitflag] = solveJumpLP(Jeq, Lxp, g, dpe, cin, y, step, z, JAbc)
% test LP with GUROBI and CLPEX
% INCLUDE ACTIVE BOUND CONSTRAINT HERE !

[numY, numX] = size(Jeq);

% active bound constraints
JActive  = [Jeq;JAbc];
numAct   = size(JAbc,1);

lb       = -Inf*ones(numY+numAct,1);
ub       = Inf*ones(numY+numAct,1);
f        = dpe*step;
f        = [f;zeros(numAct,1)];
%A        = [Jeq';-Jeq'];
A        = [JActive';-JActive'];
b        = [abs(z) - g - y.lam_x; abs(z) + g + y.lam_x];

% option   = optimoptions('linprog','Algorithm','dual-simplex','Display','off');
% option = cplexoptimset;
% option.Display = 'on';
% 
%[lpSol,~,exitflag] = linprog(f, A, b, [], [], lb, ub, y.lam_g, option );
%[lpSol,~,exitflag] = linprog(f, A, b, [], [], lb, ub, [], option );

%option  = cplexoptimset('Display', 'on', 'Algorithm', 'dual');
% reference for CLPLEX option: http://www.pserc.cornell.edu/matpower/docs/ref/matpower5.0/cplex_options.html
%option          = cplexoptimset('cplex'); % THIS CAUSES CRASH DUMP !
option          = cplexoptimset;
%option.Display  = 'iter';
option.Display  = 'none';
option.lpmethod = 1; %primal simplex
%option.lpmethod = 2; %dual simplex
%option.lpmethod = 3;  %network simplex
%option.advance  = 1;
%option.read.scale = 1;
%option.simplex.pgradient = -1; %OK - 230 sec. 116354 iterations
%option.simplex.pgradient = 3; %OK - 211 sec. 60198 iterations
option.simplex.tolerances.optimality = 1e-1; %default value is 1e-6

startlp = tic;
[lpSol,~,exitflag] = cplexlp(f, A, b, [], [], lb, ub, y.lam_g, option );

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
