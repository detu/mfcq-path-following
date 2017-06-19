function [x_init, y_init, qp_exit, delta_t, success] = solvePredictorCorrector(prob, x_init, y_init, step, lb, ub, N, x0, t, delta_t, p_0, p_t, oldEta)
%SOLVEPREDICTORCORRECTOR Summary of this function goes here
% 
% [OUTPUTARGS] = SOLVEPREDICTORCORRECTOR(INPUTARGS) Explain usage here
% Solving path-following using the standard predictor-corrector method
%
%
% Examples: 
% 
% Provide sample usage code here
% 
% See also: List related files here

% $Author: suwartad $	$Date: 2017/06/18 17:37:28 $	$Revision: 0.1 $
% Copyright: Process Control Group - NTNU Trondheim 2017

global flagDt;
flagDt = 0;
etamax  = oldEta;
success = 0;
%if oldEta > 1
if oldEta > 5
    % decrease deltaT
    delta_t = 0.6*delta_t;
    qp_exit = 0;
    % debug 
    fprintf("keyboard debug \n");
    keyboard;
else
    
    % Solve Predictor-Corrector QP
    %[deltaXp,deltaYp, qp_exit] = solveQPPredictorCorrector(H, Jeq, g, cin, Hc, step, dpe, lb, ub, deltaXc);
    [deltaXpc, ~, qp_exit, lamda, qp_run] = solveQPPredictorCorrector(prob, p_0, x_init, y_init, step, lb, ub, N, x0);
    
    % Calculate new Eta 
    xCurrent       = x_init + deltaXpc;
    yCurrent.lam_g = lamda.lam_g;
    yCurrent.lam_x = lamda.lam_x;
    flagDt         = 0;
    [~,g,~,~,cin,~,~,Jeq,~,~,~] = prob.obj(xCurrent,yCurrent,p_t, N);
    [newEta, z] = computeEta(Jeq, g, yCurrent, cin);
    
    % Checking condition (5.1)
    if newEta <= max(5.0,etamax)   % Parameters.etaMax = 5.0
        % Update primal and dual variables 
        
        success = 1;
        x_init  = xCurrent;
        y_init  = yCurrent;
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SOLVE QP program
% solution of QP program: [y] (directional derivative)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [y, qp_val, qp_exit, lamda, elapsedqp] = solveQPPredictorCorrector(prob, p, x_init, y_init, step, lb, ub, N, x0)
    

    % obtain derivatives information
    [~,g,H,Lxp,cin,~,~,Jeq,dpe,~,~] = prob.obj(x_init,y_init,p, N);    

    % QP setup
    A   = [];
    b   = [];
    f   = Lxp * step + g;
    
    % Only Equality Constraint
    ceq  = cin;
    Aeq  = Jeq;
    beq  = dpe*step + ceq;   %OK

       
%     % CHECK LAGRANGE MULTIPLIERS FROM BOUND CONSTRAINTS 
%     lmC = abs(y_init.lam_x);
%     bAc = find(lmC >= 1e-3);
%     
%     % build equality constraint from active bound constraints
%     numBaC = size(bAc,1);
%     for i=1:numBaC        
%         % put strongly active constraint on boundary
%         indB         = bAc(i);
%         
%         ub(indB)     = 0;        % keep upper bound on boundary
%         lb(indB)     = 0; 
%         
%     end

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Finally solve QP problem
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    Prob   = qpAssign(H, f, Aeq, beq, beq, lb, ub);
    Prob.optParam.eps_x = 1e-7;
    Prob.optParam.fTol  = 1e-7;
    Prob.optParam.xTol  = 1e-7;
    %Prob.PriLevOpt = 5;
    %Prob.PriLev = 1;
    startqp  = tic;
    Result = tomRun('qp-minos', Prob, 1);

    elapsedqp = toc(startqp);
    fprintf('QP solver runtime: %f\n',elapsedqp);
    qp_exit = Result.ExitFlag;
    if qp_exit == 0
        y       = Result.x_k;
        qp_val  = Result.f_k;
    else
        keyboard;
    end

    numX        = size(x_init,1);
    lamda.lam_x = Result.v_k(1:numX);
    lamda.lam_g = -Result.v_k(numX+1:end);
    fprintf('QP return: %d\n', qp_exit);
    
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