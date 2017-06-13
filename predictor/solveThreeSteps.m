function [x_init, y_init, qp_exit, delta_t, success] = solveThreeSteps(prob, x_init, y_init, step, lb, ub, N, x0, t, delta_t, p_0, p_t, oldEta)
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
if ((mpciter == 1) && (t == 0))
    flagDt = 1;
else
    flagDt = 0;
    %load Hc.mat;
end

if oldEta > 1
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
    [deltaXc,deltaYplus]      = solveCorrectStep(H, Jeq, g, cin, y_init);
    
    % Second: Predictor
    [deltaXp,deltaYp, qp_exit] = solveQPPredict(H, Jeq, g, cin, Hc, step, dpe, lb, ub, deltaXc);
    
    
    % Calculate new Eta 
    xCurrent       = x_init + deltaXc + deltaXp;
    yCurrent.lam_g = y_init.lam_g + deltaYplus + deltaYp.lam_g;
    yCurrent.lam_x = y_init.lam_x + deltaYp.lam_x;
    flagDt         = 0;
    [~,g,~,~,cin,~,~,Jeq,~,~,~] = prob.obj(xCurrent,yCurrent,p_t, N);
    [newEta, z] = computeEta(Jeq, g, yCurrent, cin);
    
    % Checking condition (5.1)
    if newEta <= max(5.0,etamax)   % Parameters.etaMax = 5.0
        % Update primal and dual variables 
        deltaX = deltaXc + deltaXp;
        deltaY = deltaYplus + deltaYp.lam_g;
        
        % Update deltaT  
        delta_t = updateDeltaT(oldEta, newEta, delta_t);
        
        % Since the dynamics are equality constraints, all constraints are
        % active. NOTE: Bound constraints will be treated later! 
        
        % Third: Jump Step
        [lpSol, exitLP]   = solveJumpLP(Jeq, Lxp, g, dpe, cin, yCurrent, step, z);
        if exitLP >= 0
            y_init.lam_g = lpSol;
            y_init.lam_x = yCurrent.lam_x; % Is it Correct? 
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

function [dXp, dYp, elapsedqp] = solveQPPredict(H, Jeq, g, cin, Hc, step, dpe, lb, ub, deltaXc)
% QP consists of equality and bound constraints 

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
% option = cplexoptimset;
% option.Display = 'on';
% 
%[lpSol,~,exitflag] = linprog(f, A, b, [], [], lb, ub, y.lam_g, option );
%[lpSol,~,exitflag] = linprog(f, A, b, [], [], lb, ub, [], option );

%option  = cplexoptimset('Display', 'on', 'Algorithm', 'dual');
% reference for CLPLEX option: http://www.pserc.cornell.edu/matpower/docs/ref/matpower5.0/cplex_options.html
option          = cplexoptimset('cplex');
%option.Display  = 'iter';
option.Display  = 'none';
option.lpmethod = 1;
option.advance  = 1;
%option.simplex.pgradient = -1; %OK - 230 sec. 116354 iterations
option.simplex.pgradient = 3; %OK - 211 sec. 60198 iterations
%option.simplex.limits.iterations = 30100;
%option.simplex.limits.iterations = 32000;

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
