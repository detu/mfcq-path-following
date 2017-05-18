function [primal,dual,deltaTc,Eta] = predictorCorrectorSLQP(caseProblem, paramInit, paramFinal, primalInit, dualInit, vargargin)
%PREDICTORCORRECTORSLQP Summary of this function goes here
% 
% Implementation of "Algorithm for Nonlinear Programming Based on Piecewise
% Linear Model" by R. Byrd, J. Nocedal, R. Waltz, and Y. Wu
%
% [OUTPUTARGS] = PREDICTORCORRECTORSLQP(INPUTARGS) Explain usage here
% 
% Examples: 
% 
% Provide sample usage code here
% 
% See also: List related files here

% $Author: suwartad $	$Date: 2017/04/24 15:20:03 $	$Revision: 0.1 $
% Copyright: Process Control Group - NTNU Trondheim 2017
%% NOTE:
% - First attempt: TRY WITH FIX STEP !

%% setup initial parameters
clear problem;
sym param;
%isComplete      = 0;
p               = paramInit;
Parameters      = setParameters();
problemName     = @(param)caseProblem(param);
problem         = problemName(p);
m               = problem.m;
me              = problem.me;
equalityIndex   = 1:me;
inequalityIndex = me+1:m;
[objectiveFunctionValue,gradientObjective] = problem.obj(primalInit,p);
[constraint,jacobianConstraint]            = problem.cons(primalInit,p);

variables.primal = primalInit;
variables.dual   = dualInit;
variables.param  = 0;
% constraints
constraints.value    = constraint;
% constraints.eqInd    = estimateActive;
% constraints.inInd    = estimateInactive;
constraints.jacobian = jacobianConstraint;
% constraints.active   = activeSet;
% objective function
objective.value      = objectiveFunctionValue;
objective.gradient   = gradientObjective;
% initialization
deltaT  = Parameters.initDeltaT;
t       = variables.param;
% p       = 0;
deltaP  = paramFinal - paramInit;
numIter = 0;  % iteration number

% initiate master dan LP trust regions
delta   = 1;
deltaLP = 0.1; % CHANGE !

while t < 1
    
    numIter = numIter + 1;
    %% compute Eta from previous step
    clear p0 problem;
    p0      = (1 - t)*paramInit + t*paramFinal;
    problem = problemName(p0);
    [objectiveFunctionValue,gradientObjective] = problem.obj(primalInit,p0);
    [constraint,jacobianConstraint]            = problem.cons(primalInit,p0);
    %% update parameter
    clear p;
    p  = (1 - t - deltaT)*paramInit + (t + deltaT)*paramFinal;
    
    % solve estimation working set by solving LP
    %stepLp  = estimateWorkingSet(gradientObjective, jacobianConstraint);
    [stepLp, penaltyParam, activeSet]  = solveLP(problem, variables, constraints, objective, p0);
    
    % compute Cauchy step and estimate Lagrange multipliers
    [xCauchy, dual] = computeCauchyPointAndDual();
    
    % compute Newton step using QP 
    xEQP    = computeEQPPoint();
    
    % compute trial point
    xTrial = computeTrialPoint();
    
    % trust-region update
    primal = trustRegionUpdate(xTrial);
    
    % Update master trust-region 
    
    % Update LP trust-region
end

end

function [stepLP, v, activeSet] = solveLP(problem, variables, constraints, objective, p)
n       = problem.n;
m       = problem.m;
me      = problem.me;
eqs     = 1:me;
iqs     = m - eqs;
% H       = problem.hess(variables.primal,variables.dual,p);
%H       = problem.hess(variables.primal,variables.dual);
J       = constraints.jacobian;
g       = objective.gradient;
% y       = variables.dual;
c       = constraints.value;
% epsA    = constraints.eqInd;

%% TO BE EDITED FROM HERE !
deltaLP = 0.8/sqrt(n);
v   = [repmat(ones(me,1),2,1) ; ones(iqs,1)];
f   = [g; v];
Aeq = [J(1:me,:) -ones(1,me) ones(1,me) zeros(1,iqs)];
beq = -c(1:me,:);
A   = -[J(me+1:end,:) zeros(iqs,me) zeros(iqs,me) diag(ones(iqs,1))];
b   = -c(me+1:end,:);
lb  = zeros(n+2*me+iqs,1);
ub  = [deltaLP*ones(n,1);inf*ones(2*me+iqs,1)]; 
%ub  = inf*ones(n+2*me+iqs,1); 

option   = optimoptions('linprog','Algorithm','dual-simplex','Display','off');
[stepLP,~,exitflag] = linprog(f, A, b, Aeq, beq, lb, ub, [], option );

% compute feasibility measure 
equality   = c(1:me) + J(1:me,:)*stepLP(1:n,:);
inequality = -c(me+1:end,:) - J(me+1:end,:)*stepLP(1:n,:);
inequality = max(0,inequality);
feasibilityMeasure = (1/m) *(sum(equality) + sum(inequality));

% update penalty parameter 
tol = 1e-8;
if feasibilityMeasure > tol
    f = [0*g; v];
    [stepLP,~,exitflag] = linprog(f, A, b, Aeq, beq, lb, ub, [], option );
    
    % compute feasibility measure
    equality   = c(1:me) + J(1:me,:)*stepLP(1:n,:);
    inequality = -c(me+1:end,:) - J(me+1:end,:)*stepLP(1:n,:);
    inequality = max(0,inequality);
    feasibilityMeasureInf = (1/m) *(sum(equality) + sum(inequality));
    
    if feasibilityMeasureInf < tol
        
        while (( feasibilityMeasure > tol ) && (v < 1e20) )
            v = 10*v;
            f = [g; v];
            [stepLP,~,exitflag] = linprog(f, A, b, Aeq, beq, lb, ub, [], option );
            
            % compute feasibility measure
            equality   = c(1:me) + J(1:me,:)*stepLP(1:n,:);
            inequality = -c(me+1:end,:) - J(me+1:end,:)*stepLP(1:n,:);
            inequality = max(0,inequality);
            feasibilityMeasure = (1/m) *(sum(equality) + sum(inequality));
        end
        
    elseif ( ( feasibilityMeasure - feasibilityMeasureInf) >= tol )
        temp = feasibilityMeasure;
        while (( (feasibilityMeasure - temp) < 0.1 * (feasibilityMeasure - feasibilityMeasureInf)) && (v < 1e20) )
            v = 10*v;
            f = [g; v];
            [stepLP,~,exitflag] = linprog(f, A, b, Aeq, beq, lb, ub, [], option );
            
            % compute feasibility measure
            equality   = c(1:me) + J(1:me,:)*stepLP(1:n,:);
            inequality = -c(me+1:end,:) - J(me+1:end,:)*stepLP(1:n,:);
            inequality = max(0,inequality);
            temp = (1/m) *(sum(equality) + sum(inequality));
        end
        
    end
    
end

% COMPUTE ACTIVE SET

end

function [xC, dual] = computeCauchyPointAndDual()

end

function xQ = computeEQPPoint()

end

function xTrial = computeTrialPoint()

end

function primal = trustRegionUpdate()

end

function Parameters = setParameters()
Parameters.gamma        = 0.7;
Parameters.alpha        = 0.6;
%Parameters.alpha        = 0.9;
Parameters.maxIteration = 50000;
%Parameters.initDeltaT   = 0.05;
Parameters.initDeltaT   = 0.1;
Parameters.maxFailure   = 50;
%Parameters.maxFailure   = 40;
Parameters.optValue     = 1e-6;
Parameters.etaMin       = 1e-6;
%Parameters.etaMax       = 1e-2;
Parameters.etaMax       = 1e-4;

end