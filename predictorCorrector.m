function [primal,dual] = predictorCorrector(caseProblem, paramInit, paramFinal, primalInit, dualInit, vargargin)
%PREDICTORCORRECTOR Summary of this function goes here
%
% [OUTPUTARGS] = PREDICTORCORRECTOR(INPUTARGS) Explain usage here
%
% Implementation of Predictor-Corrector QP for MFCQ case
%
% Examples:
%
% Provide sample usage code here
%
% See also: List related files here

% $Author: suwartad $    $Date: 2017/04/04 17:44:45 $    $Revision: 0.1 $
% Copyright: Process Control Group - NTNU Trondheim 2017

%% setup parameters
clear problem;
sym param;
p               = paramInit;
Parameters      = setParameters();
problemName     = @(param)caseProblem(param);
problem         = problemName(p);
m               = problem.m;
me              = problem.me;
%eqs             = 1:me;
% ineqs           = me+1:m;
% n               = problem.n ;
% numEquality     = problem.m;
% numInequality   = problem.me;
equalityIndex   = 1:me;
inequalityIndex = me+1:m;
numFail         = 0;             % set counter for fail attempt
%% initial estimate of Active Set
[objectiveFunctionValue,gradientObjective] = problem.obj(primalInit,p);
[constraint,jacobianConstraint]            = problem.cons(primalInit,p);
[currentEta,~]                             = calculateEta(problem, dualInit, gradientObjective, constraint, jacobianConstraint);
[activeSet, estimateActive,estimateInactive] = estimateActiveSet(dualInit, constraint, equalityIndex, inequalityIndex, currentEta, Parameters);

etaMax = Parameters.etaMax;
% make objects in compact form
% variables (primal, dual, and parameter t) -> MAKE AS A FUNCTION !!!
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
primal  = [];
dual    = [];
deltaTc = [];
Eta     = [];
tData   = [];

% main loop
while t < 1
    numIter = numIter + 1;
    %% update parameter
    clear p problem;
    p  = (1 - t - deltaT)*paramInit + (t + deltaT)*paramFinal;
    problem = problemName(p);
    [objectiveFunctionValue,gradientObjective] = problem.obj(primalInit,p);
    [constraint,jacobianConstraint]            = problem.cons(primalInit,p);
    % constraints
    constraints.value    = constraint;
    constraints.eqInd    = estimateActive;
    constraints.inInd    = estimateInactive;
    constraints.jacobian = jacobianConstraint;
    constraints.active   = activeSet;
    % objective function
    objective.value      = objectiveFunctionValue;
    objective.gradient   = gradientObjective;
    [currentEta,~]                             = calculateEta(problem, dualInit, gradientObjective, constraint, jacobianConstraint);
    
    if currentEta >= 1
        deltaT  = deltaT * Parameters.alpha;
        p       = (1 - t - deltaT)*paramInit + (t + deltaT)*paramFinal;
        clear problem;
        problem = problemName(p);
        [objectiveFunctionValue,gradientObjective] = problem.obj(primalInit,p);
        [constraint,jacobianConstraint]            = problem.cons(primalInit,p);
        % constraints
        constraints.value    = constraint;
        constraints.eqInd    = estimateActive;
        constraints.inInd    = estimateInactive;
        constraints.jacobian = jacobianConstraint;
        constraints.active   = activeSet;
        % objective function
        objective.value      = objectiveFunctionValue;
        objective.gradient   = gradientObjective;
    end
    
    %% SOLVE PC QP
    [deltaX, deltaY, exitQP] = solveQPPredictorCorrector(problem, objective, variables, constraints, p);
    clear problem;
    problem = problemName(p);
    [~,gradientObjective]           = problem.obj(primalInit+deltaX,p);
    [constraint,jacobianConstraint] = problem.cons(primalInit+deltaX,p);
%     % constraints
%     constraints.value    = constraint;
%     constraints.eqInd    = estimateActive;
%     constraints.inInd    = estimateInactive;
%     constraints.jacobian = jacobianConstraint;
%     constraints.active   = activeSet;
%     % objective function
%     objective.value      = objectiveFunctionValue;
%     objective.gradient   = gradientObjective;
    [nextEta,Lag]  = calculateEta(problem, deltaY, gradientObjective, constraint, jacobianConstraint);
    
    %while( nextEta >= max(currentEta,Parameters.etaMax) )
    %while( nextEta > Parameters.etaMax )
    while( nextEta > Parameters.etaMax ) || ( nextEta > etaMax)
        % decrease delta_t
        deltaT  = deltaT * Parameters.alpha;
        p       = (1 - t - deltaT)*paramInit + (t + deltaT)*paramFinal;
        numFail = numFail + 1;
        %clear deltaX deltaY problem;
        %problem = problemName(p0);
        %clear problem;
        %problem = problemName(p);
        [objectiveFunctionValue,gradientObjective] = problem.obj(primalInit,p);
        [constraint,jacobianConstraint]            = problem.cons(primalInit,p);
        % constraints
        constraints.value    = constraint;
        constraints.eqInd    = estimateActive;
        constraints.inInd    = estimateInactive;
        constraints.jacobian = jacobianConstraint;
        constraints.active   = activeSet;
        % objective function
        objective.value      = objectiveFunctionValue;
        objective.gradient   = gradientObjective;
        
        % solve QPPredict again  --> SOLVE PREDICTOR CORRECTOR QP HERE !
        %[deltaXp,deltaYp,exitQP] = solveQPPredict(problem, variables, constraints, deltaXc, deltaT, deltaP);
        [deltaX, deltaY, exitQP] = solveQPPredictorCorrector(problem, objective, variables, constraints, p);
        %deltaX  = deltaXc + deltaXp;
        %deltaY  = deltaYplus + deltaYp;
        %clear problem;
        %problem = problemName(p);
        [~,gradientObjective]           = problem.obj(primalInit+deltaX,p);
        [constraint,jacobianConstraint] = problem.cons(primalInit+deltaX,p);
%         % constraints
%         constraints.value    = constraint;
%         constraints.eqInd    = estimateActive;
%         constraints.inInd    = estimateInactive;
%         constraints.jacobian = jacobianConstraint;
%         constraints.active   = activeSet;
%         % objective function
%         objective.value      = objectiveFunctionValue;
%         objective.gradient   = gradientObjective;
        
        %[nextEta,Lag]  = calculateEta(problem, dualInit+deltaY, gradientObjective, constraint, jacobianConstraint);
        [nextEta,Lag]  = calculateEta(problem, deltaY, gradientObjective, constraint, jacobianConstraint);
    end
    
    % check good eta condition, if satisfied increase delta_t
    if (nextEta < currentEta^(1+Parameters.gamma))
    %if (nextEta < currentEta)
        deltaT = min(1-t,deltaT/Parameters.alpha);
    else
        deltaT = min(1-t,deltaT);
    end

    % heuristic to prevent small Eta
    if nextEta < 1e-8
        deltaT = min(1-t,deltaT/Parameters.alpha);
    end

    % compute estimate of active set A
    currentEta = nextEta;
    if currentEta > etaMax
        etaMax = currentEta;
    end
    dualInit   = dualInit+deltaY;
    primalInit = primalInit+deltaX;
    [activeSet, estimateActive,estimateInactive] = estimateActiveSet(dualInit, constraint, equalityIndex, inequalityIndex, currentEta, Parameters);

    %% SOLVE JumpLP
    variables.primal = primalInit;
    variables.dual   = dualInit;
    % constraints
    constraints.value    = constraint;
    constraints.eqInd    = estimateActive;
    constraints.inInd    = estimateInactive;
    constraints.jacobian = jacobianConstraint;
    constraints.active   = activeSet;
    % objective function
    objective.value      = objectiveFunctionValue;
    objective.gradient   = gradientObjective;
    [dualLP,exitLP]      = solveLP(problem, objective, variables, constraints, deltaT, deltaP, Lag);

    linobj = problem.dcdp(activeSet,:)*deltaT*deltaP;
    % define A+ set
    if exitLP > 0 && (linobj'*dualLP < linobj'*dualInit(activeSet)-sqrt(eps))
    %if exitLP > 0
        dualInit(activeSet) = dualLP;
        [activeSet, estimateActive,estimateInactive] = estimateActiveSet(dualInit, constraint, equalityIndex, inequalityIndex, currentEta, Parameters);
        %variables.primal = primalInit;
        variables.dual   = dualInit;
        % constraints
        %constraints.value    = constraint;
        constraints.eqInd    = estimateActive;
        constraints.inInd    = estimateInactive;
        %constraints.jacobian = jacobianConstraint;
        constraints.active   = activeSet;
    end

    % update t
    t          = t + deltaT;
    currentEta = nextEta;
    primal     = [primal;primalInit'];
    dual       = [dual;dualInit'];
    deltaTc    = [deltaTc;deltaT];
    Eta        = [Eta;currentEta];
    tData      = [tData;t];
    % print out every iteration ...
    fprintf('Iter = %6g    Eta = %5.2e    t =  %3.1e    deltaT = %3.1e    QPflag=%d \n', numIter, currentEta, t, deltaT, exitQP);
        
end
%% PLOT PRIMAL, DUAL, ETA, and deltaT !!!
keyboard;

end

function Parameters = setParameters()
Parameters.gamma        = 0.7;
Parameters.alpha        = 0.6;
%Parameters.alpha        = 0.9;
Parameters.maxIteration = 50000;
%Parameters.initDeltaT   = 0.01;
%Parameters.initDeltaT   = 0.001;
Parameters.initDeltaT   = 0.1; % OK
Parameters.maxFailure   = 50;
%Parameters.maxFailure   = 40;
Parameters.optValue     = 1e-6;
Parameters.etaMin       = 1e-6;
%Parameters.etaMax       = 1e-2;
Parameters.etaMax       = 1e-4; % OK
%Parameters.etaMax       = 1e-6;
end

function [eta,Lagrangian] = calculateEta(problem, dual, grad, cons, Jac)
% compute Lagrangian
Lagrangian = grad - Jac'*dual;

% equality constraints
equalityConstraints   = cons(1:problem.me);
inequalityConstraints = cons(problem.me+1:end);

% bigger value between inequality constraints and dual variables
numInequalityConstraints = size(inequalityConstraints,1);
minInequalityAndDual     = zeros(numInequalityConstraints,1);
correspondingDual        = dual(problem.me+1:end);
for i=1:numInequalityConstraints
    minInequalityAndDual(i,1) = min(inequalityConstraints(i),correspondingDual(i));
end

% construct vector of Eta
vectorEta = [Lagrangian;equalityConstraints;minInequalityAndDual];
% Eta is infinite norm of that vector
eta       = norm(vectorEta,inf);
%eta       = 0.1*norm(vectorEta,inf);

end

function [activeSet, estimateActive, estimateInactive]  = estimateActiveSet(dualInit, constraint, equalityIndex, inequalityIndex, Eta, Parameters)
estimateActive   = find(abs(dualInit) > Eta^Parameters.gamma);
estimateActive   = setdiff(estimateActive,equalityIndex);
estimateInactive = setdiff(inequalityIndex ,estimateActive);
activeConstraint = union(equalityIndex,find(constraint < Eta^(Parameters.gamma)));
activeSet        = union(estimateActive,activeConstraint);
end

function [deltaXp,deltaYp, exitflag] = solveQPPredictorCorrector(problem,  objective, variables, constraints, p)
n       = problem.n;
m       = problem.m;
me      = problem.me;
eqs     = 1:me;
H       = problem.hess(variables.primal,variables.dual,p);
%H       = problem.hess(variables.primal,variables.dual);
% dcdp    = problem.dcdp;
J       = constraints.jacobian;
g       = objective.gradient;
% y       = variables.dual;
c       = constraints.value;
epsA    = constraints.eqInd;
epsF    = constraints.inInd;


% [deltax,~,exitflag,~,lambda] = quadprog(H,g,-J(epsF,:),c(epsF),[J(eqs,:);J(epsA,:)],[-c(eqs,:);-c(epsA,:)],[],[], x,options);

A   = -J(epsF,:);
b   = c(epsF);
Aeq = [J(eqs,:);J(epsA,:)];
beq = [-c(eqs,:);-c(epsA,:)];
option = optimset('Display','off','Algorithm','active-set');
[deltaXp,~,exitflag,~,lambda] = quadprog(H,g,A,b,Aeq,beq,[],[], variables.primal, option);
% deltaYp             = zeros(m,1);
% deltaYp(epsF)       = (lambda.ineqlin(1:(length(epsF))));
% deltaYp([eqs;epsA]) = -lambda.eqlin;
% NEED TO CHECK exitflag HERE !!!
if exitflag < 0
    keyboard;
    deltaXp = zeros(n,1);
    deltaYp = zeros(m,1);
else
    deltaYp             = zeros(m,1);
    deltaYp(epsF)       = (lambda.ineqlin(1:(length(epsF))));
    deltaYp([eqs;epsA]) = -lambda.eqlin;
end
end

function [y,exitflag] = solveLP(problem, objective, variables, constraints, deltaT, deltaP, z)

Act  = constraints.active;
me   = problem.me;
dcdp = problem.dcdp;
J    = constraints.jacobian;
g    = objective.gradient;

lb       = zeros(length(Act),1);
lb(1:me) = -Inf;
ub       = Inf*ones(length(Act),1);
f        = dcdp(Act,:)*deltaT*deltaP;
A        = [J(Act,:)';-J(Act,:)'];
b        = [g+abs(z);abs(z)-g];
option   = optimoptions('linprog','Algorithm','dual-simplex','Display','off');
[y,~,exitflag] = linprog(f, A, b, [],[],lb,ub,variables.dual(Act),option );
end