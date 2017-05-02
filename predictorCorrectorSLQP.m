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
constraints.eqInd    = estimateActive;
constraints.inInd    = estimateInactive;
constraints.jacobian = jacobianConstraint;
constraints.active   = activeSet;
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
    
    % 1. solve estimation working set by solving LP
    %stepLp  = estimateWorkingSet(gradientObjective, jacobianConstraint);
    stepLp  = estimateWorkingSet(problem, variables, constraints, objective, p0);
    
    % 2. compute Cauchy step and estimate Lagrange multipliers
    [xCauchy, dual] = computeCauchyPointAndDual();
    
    % 3. compute Newton step using QP 
    xEQP    = computeEQPPoint();
    
    % 4. compute trial point
    xTrial = computeTrialPoint();
    
    % 5. trust-region update
    primal = trustRegionUpdate(xTrial);
    
    % 6. Update master trust-region 
    
    % 7. Update LP trust-region
end

end

function stepLP = estimateWorkingSet(problem, variables, constraints, objective, p)
n       = problem.n;
m       = problem.m;
me      = problem.me;
eqs     = 1:me;
H       = problem.hess(variables.primal,variables.dual,p);
%H       = problem.hess(variables.primal,variables.dual);
J       = constraints.jacobian;
g       = objective.gradient;
y       = variables.dual;
c       = constraints.value;
epsA    = constraints.eqInd;


end

function [xC, dual] = computeCauchyPointAndDual()

end

function xQ = computeEQPPoint()

end

function xTrial = computeTrialPoint()

end

function primal = trustRegionUpdate()

end