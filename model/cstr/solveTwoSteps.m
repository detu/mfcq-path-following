function [x_init, y_init, qp_exit, delta_t, success, etaData, numActiveBound, activeBoundTol] = solveTwoSteps(prob, x_init, y_init, step, lb, ub, N, x0, t, delta_t, p_0, p_t, oldEta, ub_init, gOld)
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
% Implement 2 steps: Predictor-corrector QP and LP.
% Remove second order (HESSIAN of CONSTRAINT) information
% Increase eta (KKT residual tolerance) to 0.1 (or bigger)

% Derivatis w.r.t to deltaT can be computed ONE TIME.
global flagDt mpciter;
etamax         = oldEta;
global success;
keepgoingrank  = 1;      % from Slava Code  
%keepgoingrank = 0;
prevRankDef    = 0;
count          = 0;
flagDt         = 0;
numX           = size(x_init,1);
numY           = size(y_init.lam_g,1);
activeBoundTol = [];

% % %if oldEta > 1
% if oldEta > 5
% % if oldEta > 100
%     % decrease deltaT
%     delta_t = 0.6*delta_t;
%     qp_exit = 0;
%     % debug 
%     fprintf('keyboard debug \n');
%     keyboard;
%     
%     success = 0;
%     etaData        = [];
%     numActiveBound = [];
%     activeBoundTol = [];
% else
    % proceed with everything... 
    
    
    % obtain derivatives information
    [~,g,H,Lxp,cin,~,~,Jeq,dpe,~,~] = prob.obj(x_init,y_init,p_0, N);

    % First step: Predictor-Corrector QP
    [deltaXpc,deltaYpc, qp_exit, JAbc, activeBoundTol, activeIndex] = solvePredictorCorrectorQP(H, Jeq, g, gOld, cin, step, dpe, lb, ub, y_init, oldEta);
    
    % feasible QP
    if qp_exit == 0
        
        % Calculate new Eta
        xCurrent         = x_init + deltaXpc;
        yCurrent.lam_g   = deltaYpc.lam_g;
        
        % update the dual's bound constraints
        yCurrent.lam_x   = deltaYpc.lam_x;
        
        flagDt         = 0;
        [~,g,~,~,cin,~,~,Jeq,~,~,~] = prob.obj(xCurrent,yCurrent,p_t, N);
        [newEta, z] = computeEta(Jeq, g, yCurrent, cin, activeBoundTol);
        fprintf('newEta = %e \n',newEta);
        
        % Checking condition (5.1)
        %if newEta <= 1e-1
        if newEta <= 1e9
            
            % PREDICTOR-CORRECTOR QP
            deltaY = deltaYpc.lam_g;
            y_init.lam_g = deltaYpc.lam_g;
            y_init.lam_x = deltaYpc.lam_x;
            
            % update primal variable
            x_init  = xCurrent;
            
            % Update deltaT
            delta_t = updateDeltaT(oldEta, newEta, delta_t);
            fprintf('new delta_t =%f\n', delta_t);
            
            % solve LP only for MFCQ
            %numActiveBound = size(activeBoundInd,1);
            %if numY+numActiveBound > numX
                
                % Second step: Jump Step
                [lpSol, exitLP, linobj]   = solveJumpLP(Jeq, Lxp, g, dpe, cin, yCurrent, step, z, JAbc, activeIndex);
                
                if exitLP >= 0
                    %condLP1 = linobj'*lpSol;
                    condLP1 = linobj'*-lpSol;
                    condLP2 = linobj'*[deltaY;y_init.lam_x(activeIndex)]-sqrt(eps);
                    if (condLP1 < condLP2)
                    %if (condLP1 > condLP2)
                        %y_init.lam_g = lpSol(1:numY);
                        y_init.lam_g = -lpSol(1:numY);
                        y_init.lam_x = zeros(numX,1);
                        y_init.lam_x(activeIndex) = -lpSol(numY+1:end); % CHECK NEGATIF OR POSITIF ?
                        %y_init.lam_x(activeBoundInd) = lpSol(numY+1:end);
                    end
                end
            %end
            
            success = 1;
            qp_exit = 1;
            
            % Update etaMax
            if newEta > etamax
                etamax = newEta;
            end
            
            etaData        = newEta;
            numActiveBound = numel(activeIndex);
        else
            fprintf('delta_t =%f\n', delta_t);
            % decrease deltaT
            delta_t = 0.6*delta_t;
            qp_exit = 0;
            % debug
            fprintf('keyboard debug \n');
            if (delta_t <= 1e-3)
                keyboard;
            end
            success = 0;
            etaData        = [];
            numActiveBound = [];
            activeBoundTol = [];
            %flagDt = 0;
            
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



function [dXp, dYp, qp_exit, JAbc, activeBoundTol, activeIndex] = solvePredictorCorrectorQP(Hm, Jeq, g, gOld, cin, step, dpe, lb, ub, y, oldEta)

% active bound constraints
boundMultiplier    = y.lam_x;
positiveBoundMult  = abs(boundMultiplier);
activeIndex        = find(positiveBoundMult>1e-1);  % set active bound constraint.
activeBoundTol     = oldEta^(0.7);
% activeIndex        = find(positiveBoundMult>activeBoundTol);  % set active bound constraint.
%paramIndex         = find(activeIndex <= 84);       % remove the first 84 constraints (for Distillation column!)
%activeIndex(paramIndex) = [];
numActiveBoundCons = numel(activeIndex);

[numY,m] = size(Jeq);
JAbc     = zeros(numActiveBoundCons,m);
for i=1:numActiveBoundCons
    row         = activeIndex(i);
    JAbc(i,row) = 1;
end

% QP setup
%A   = [];
%b   = [];
%f   = [];
%f   = g;
f = g - gOld;

% Only Equality Constraint
ceq         = cin;
[numY,numX] = size(Jeq);

Aeq  = Jeq;
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
startqp  = tic;
Result = tomRun('qp-minos', Prob, 1);

elapsedqp = toc(startqp);
fprintf('QP solver runtime: %f\n',elapsedqp);
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

function [lpSol, exitflag, f] = solveJumpLP(Jeq, Lxp, g, dpe, cin, y, step, z, JAbc, activeIndex)


[numY, numX]       = size(Jeq);

% active bound constraints
JActive  = [Jeq;JAbc];
numAct   = size(JAbc,1);

% lb       = -Inf*ones(numY+numAct,1);
% ub       = Inf*ones(numY+numAct,1);
lb       = [];
ub       = [];
f        = dpe*step;
f        = [f;zeros(numAct,1)];
%A        = [Jeq';-Jeq'];
A        = [JActive';-JActive'];
% %b        = [abs(z) - g - y.lam_x; abs(z) + g + y.lam_x]; % WRONG!
% b        = [abs(z) - g - JActive'*([y.lam_g; y.lam_x(activeIndex)]); abs(z) + g + JActive'*([y.lam_g; y.lam_x(activeIndex)])];
% %b        = [abs(z) + g + JActive'*([y.lam_g; y.lam_x(activeIndex)]); abs(z) - g - JActive'*([y.lam_g; y.lam_x(activeIndex)])];
b        = [abs(z) - g - JActive'*([y.lam_g; y.lam_x(activeIndex)]); abs(z) + g + JActive'*([y.lam_g; y.lam_x(activeIndex)])];

% % try with different constraint
% Aeq = JActive';
% beq = - g - JActive'*([y.lam_g; y.lam_x(activeIndex)]);

% lb       = -Inf*ones(numY+numAct,1);
% ub       = Inf*ones(numY+numAct,1);
% f        = dpe*step;
% f        = [f;zeros(numAct,1)];
% %A        = [Jeq';-Jeq'];
% %b        = [abs(z) - g - Jeq'*y.lam_g; abs(z) + g + Jeq'*y.lam_g];
% A        = JActive';
% b        = - g - JActive'*([y.lam_g; y.lam_x(activeIndex)]);

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
%option.simplex.tolerances.optimality = 1e-6; %default value is 1e-6

startlp = tic;
[lpSol,~,exitflag] = cplexlp(f, A, b, [], [], lb, ub, y.lam_g, option );
%[lpSol,~,exitflag] = cplexlp(-f, A, b, [], [], lb, ub, y.lam_g, option ); %maximization problem
%[lpSol,~,exitflag] = cplexlp(f, [], [], A, b, lb, ub, y.lam_g, option );

elapsedlp = toc(startlp);
% lpSol     = -lpSol; % different sign !
fprintf('LP solver runtime: %f\n',elapsedlp);
fprintf('LP exitflag: %d\n', exitflag);
% Take out active bound constraint here !


end

function deltaT = updateDeltaT(currentEta, nextEta, t)
gamma = 0.7;
alpha = 0.6;

% check good eta condition, if satisfied increase delta_t
if (nextEta < currentEta^(1+gamma))
    %if (nextEta < currentEta)
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

