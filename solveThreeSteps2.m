function [x_init, y_init, runtime, delta_t, success, etaData, numActiveBound, activeBoundTol] = solveThreeSteps2(prob, x_init, y_init, step, lb, ub, N, x0, t, delta_t, p_0, p_t, oldEta, ub_init, gOld)
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
etamax        = oldEta;
global success;
keepgoingrank = 1;      % from Slava Code
prevRankDef   = 0;
count         = 0;
if (t == 0) && (success == 1)
    flagDt = 1;
else
    flagDt = 0;
    %load Hc.mat;
end
numX = size(x_init,1);
numY = size(y_init.lam_g,1);


% obtain derivatives information
[~,g,H,Lxp,cin,~,~,Jeq,dpe,~,Hc] = prob.obj(x_init,y_init,p_0, N);
if isempty(Hc)
    load Hc.mat;    % CSTR + Dist. Col. A
    %load HcCstr.mat; % CSTR only
end

% debug code
%ydebug = y_init;

% keepgoingrank routine from Slava's code
if t~= 0
    %y_modif = y_init;
    while keepgoingrank
        %[y_modif, keepgoingrank, prevRankDef, count] = updateDualBasedOnRank(y_modif, Jeq, prevRankDef, count, oldEta);
        [y_init, keepgoingrank, prevRankDef, count] = updateDualBasedOnRank(y_init, Jeq, prevRankDef, count, oldEta);
        if keepgoingrank == 0
            break;
        end
    end
    % First: Corrector Step
    [deltaXc,deltaYplus,activeBoundTol,correctRt]      = solveCorrectStep(H, Hc, Jeq, g, cin, y_init, x_init, ub_init, oldEta, t);
    %[deltaXc,deltaYplus,activeBoundTol,correctRt]      = solveCorrectStep(H, Hc, Jeq, g, cin, y_modif, x_init, ub_init, oldEta, t);
else
    % skip Corrector Step at t=0 (beginning of path-following)
    [~,~,activeBoundTol,correctRt]      = solveCorrectStep(H, Hc, Jeq, g, cin, y_init, x_init, ub_init, oldEta, t);
    deltaXc          = zeros(numX,1);
    deltaYplus.lam_g = zeros(numY,1);
    deltaYplus.lam_x = zeros(numX,1);
end
    

% Second: Predictor
[deltaXp,deltaYp, qp_exit, predictRt] = solveQPPredict(H, Jeq, g, gOld, cin, Hc, step, dpe, lb, ub, deltaXc);

% feasible QP
if qp_exit == 0
    
    % Calculate new Eta
    xCurrent         = x_init + deltaXc + deltaXp;
    yCurrent.lam_g   = y_init.lam_g + deltaYplus.lam_g + deltaYp.lam_g;
    
    % update the dual's bound constraints
    yCurrent.lam_x = y_init.lam_x + deltaYplus.lam_x + deltaYp.lam_x;
    
    flagDt         = 0;
    [~,g,~,~,cin,~,~,Jeq,~,~,~] = prob.obj(xCurrent,yCurrent,p_t, N);
    [newEta, z] = computeEta(Jeq, g, yCurrent, cin, activeBoundTol);
    fprintf('newEta = %e \n',newEta);
    
    % Checking condition (5.1)
    %if newEta <= max(5,etamax)   % Parameters.etaMax = 5.0
    %if newEta <= 50
    if newEta <= 1e6
        %if newEta <= 2e-3
        deltaY = deltaYplus.lam_g + deltaYp.lam_g;
        y_init.lam_g = yCurrent.lam_g;
        y_init.lam_x = yCurrent.lam_x;
        
        % update primal variable
        x_init  = xCurrent;
        
        % Update deltaT
        delta_t = updateDeltaT(oldEta, newEta, delta_t);
        fprintf('new delta_t =%f\n', delta_t);
        
        % Third: Jump Step
        % compute new Y / dual variable here !
        [lpSol, exitLP, linobj, activeBoundInd, jumpRt]   = solveJumpLP(Jeq, g, dpe, yCurrent, step, z);
        
        if exitLP >= 0
            %condLP1 = linobj'*lpSol;
            condLP1 = linobj'*-lpSol;
            condLP2 = linobj'*[deltaY;y_init.lam_x(activeBoundInd)]-sqrt(eps);
            if (condLP1 < condLP2)
                %if (condLP1 > condLP2)
                %y_init.lam_g = lpSol(1:numY);
                y_init.lam_g = -lpSol(1:numY);
                y_init.lam_x = zeros(numX,1);
                %y_init.lam_x(activeBoundInd) = -lpSol(numY+1:end); % CHECK NEGATIF OR POSITIF ?
                %y_init.lam_x(activeBoundInd) = lpSol(numY+1:end);
                y_init.lam_x(activeBoundInd) = abs(lpSol(numY+1:end));
                % delete negative multipliers
                %negLamX                   = find(y_init.lam_x<0);
                %dYplus.lam_x(negLamX)     = 0;
            end
        end
        
        success = 1;
        qp_exit = 1;
        runtime.correctRt = correctRt;
        runtime.predictRt = predictRt;
        runtime.jumpRt    = jumpRt;
        %runtime.jumpRt    = 0;    % without jumpLP
        
        % Update etaMax
        if newEta > etamax
            etamax = newEta;
        end
        
        etaData        = newEta;
        numActiveBound = numel(activeBoundInd);
    else
        fprintf('delta_t =%f\n', delta_t);
        % decrease deltaT
        delta_t = 0.6*delta_t;
        qp_exit = 0;
        % debug
        fprintf('keyboard debug \n');
        keyboard;
        if (delta_t <= 1e-3)
            keyboard;
        end
        success = 0;
        etaData        = [];
        numActiveBound = [];
        activeBoundTol = [];
        %flagDt = 0;
        runtime = 0;
        
    end
else % infeasible QP
    fprintf('delta_t =%f\n', delta_t);
    % decrease deltaT
    delta_t = 0.6*delta_t;
    qp_exit = 0;
    % debug
    fprintf('keyboard debug \n');
    if (delta_t <= 1e-3)
        keyboard;
    end
    %keyboard;
    success = 0;
    etaData        = [];
    numActiveBound = [];
    activeBoundTol = [];
end

end

function [dXc,dYplus,activeBoundTol,correctRt] = solveCorrectStep(H, Hc, Jeq, g, c, y, x, ub, oldEta, t)

% active bound constraints
boundMultiplier    = y.lam_x;
positiveBoundMult  = abs(boundMultiplier);
activeIndex        = find(positiveBoundMult>1e-1);  % set active bound constraint.
activeBoundTol     = oldEta^(0.7);
numActiveBoundCons = numel(activeIndex);

if t == 0
    dXc       = 0;
    dYplus    = 0;
    correctRt = 0;
else
    [numY,m] = size(Jeq);
    JAbc     = zeros(numActiveBoundCons,m);
    for i=1:numActiveBoundCons
        row         = activeIndex(i);
        JAbc(i,row) = 1;
    end
    
    Jm = [Jeq;JAbc];
    n  = size(Jm,1);
    Hm = H;
    
    lhs     = [Hm                   Jm'; ...
        Jm                   zeros(n,n)];
    rhs     = -[g + Jeq'*y.lam_g + JAbc'*y.lam_x(activeIndex); [c;zeros(numActiveBoundCons,1)]];
    
    % solve a system of linear equations
    timeLE          = tic;
    solCorrectStep  = lhs\rhs;
    correctRt       = toc(timeLE);
    fprintf('Corrector Step -  solver runtime: %f\n', correctRt);
    % collect results
    dXc                       = solCorrectStep(1:m);
    dYplus.lam_g              = solCorrectStep(m+1:m+numY);
    dYplus.lam_x              = zeros(m,1);
    dYplus.lam_x(activeIndex) = solCorrectStep(m+numY+1:end);
    % delete negative multipliers
    negLamX                   = find(dYplus.lam_x<0);
    dYplus.lam_x(negLamX)     = 0;
end
end

function [dXp, dYp, qp_exit, elapsedqp] = solveQPPredict(Hm, Jeq, g, gOld, cin, Hc, step, dpe, lb, ub, deltaXc)


% QP setup
%A   = [];
%b   = [];
%f   = [];
%f   = g;
f = g - gOld;

% Only Equality Constraint
ceq         = cin;
[numY,numX] = size(Jeq);
Hcx         = zeros(numY,numX);
for i=1:numY
    Hcx(i,:) = (Hc{:,:,i}*deltaXc)';
end
Aeq  = Jeq + Hcx;
% Aeq  = Jeq;
beq  = dpe*step + ceq;   %OK

% % build equality constraint from active bound constraints
% numBaC = size(activeBoundInd,1);
% for i=1:numBaC
%     % put strongly active constraint on boundary
%     indB         = activeBoundInd(i);
%     %ub(indB)     = 0;        % keep upper bound on boundary
%     %lb(indB)     = 0;
%     %lb(indB)     = -5e-4;   %OK
%     lb(indB)     = -1e-3;
% end

% TOMLAB setup
% CHANGE H here ! H should include the Hessian of the constraints ! (same
% as in the Corrector Step)
Prob   = qpAssign(Hm, f, Aeq, beq, beq, lb, ub);

Prob.optParam.eps_x = 1e-7;
Prob.optParam.fTol  = 1e-7;
Prob.optParam.xTol  = 1e-7;
%startqp  = tic;
Result = tomRun('qp-minos', Prob, 1);

%elapsedqp = toc(startqp);
elapsedqp = Result.REALtime;
fprintf('QPPredict Step - solver runtime: %f\n',elapsedqp);
qp_exit = Result.ExitFlag;
if qp_exit == 0
    dXp       = Result.x_k;
    %qp_val  = Result.f_k;
    dYp.lam_x = -Result.v_k(1:numX);
    dYp.lam_g = -Result.v_k(numX+1:end);
    
else
    %keyboard;
    dXp       = zeros(numX,1);
    dYp.lam_x = zeros(numX,1);
    dYp.lam_g = zeros(numY,1);
end

end

function [lpSol, exitflag, f, activeIndex, elapsedlp] = solveJumpLP(Jeq, g, dpe, y, step, z)

% active bound constraints
boundMultiplier    = y.lam_x;
positiveBoundMult  = abs(boundMultiplier);
activeIndex        = find(positiveBoundMult>1e-1);  % set active bound constraint.
%activeIndex        = find(positiveBoundMult>1e-3); 
numActiveBoundCons = numel(activeIndex);

[~,m] = size(Jeq);
JAbc     = zeros(numActiveBoundCons,m);
for i=1:numActiveBoundCons
    row         = activeIndex(i);
    JAbc(i,row) = 1;
end

JActive  = [Jeq;JAbc];
numAct   = size(JAbc,1);
lb       = [];
ub       = [];
f        = dpe*step;
f        = [f;zeros(numAct,1)];
A        = [JActive';-JActive'];
b        = [abs(z) - g - JActive'*([y.lam_g; y.lam_x(activeIndex)]); abs(z) + g + JActive'*([y.lam_g; y.lam_x(activeIndex)])];
%b        = [abs(z) - g ; abs(z) + g ];


%option  = cplexoptimset('Display', 'on', 'Algorithm', 'dual');
% reference for CLPLEX option: http://www.pserc.cornell.edu/matpower/docs/ref/matpower5.0/cplex_options.html
%option          = cplexoptimset('cplex'); % THIS CAUSES CRASH DUMP !
option          = cplexoptimset;
%option.Display  = 'iter';
option.Display  = 'none';
option.lpmethod = 1; %primal simplex

startlp = tic;
%[lpSol, ~, exitflag, output] = cplexlp(f, A, b, [], [], lb, ub, y.lam_g, option );
[lpSol, ~, exitflag, output] = cplexlp(f, A, b, [], [], lb, ub, [], option );
elapsedlp = toc(startlp);
% lpSol     = -lpSol; % different sign !
fprintf('JumpLP Step - solver runtime: %f\n',elapsedlp);
fprintf('LP exitflag: %d\n', exitflag);
% Take out active bound constraint here !


end

function deltaT = updateDeltaT(currentEta, nextEta, t)
gamma = 0.7;
alpha = 0.6;

% check good eta condition, if satisfied increase delta_t
if (nextEta < currentEta^(1+gamma))
    %if (nextEta < currentEta)
    %deltaT = min(1-t,t/alpha);
    if t == 1
        deltaT = 1;
    else
        deltaT = min(1-t,t/alpha);
    end
else
    if t == 1
        deltaT = 1;
    else
        deltaT = min(1-t,t);
    end
end

end
