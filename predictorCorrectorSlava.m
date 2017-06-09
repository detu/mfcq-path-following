function [primal,dual,deltaTc,Eta] = predictorCorrectorSlava(caseProblem, paramInit, paramFinal, primalInit, dualInit, vargargin)
%PREDICTORCORRECTORSLAVA Summary of this function goes here
% 
% [OUTPUTARGS] = PREDICTORCORRECTORSLAVA(INPUTARGS) Explain usage here
% 
% Examples: 
% 
% Provide sample usage code here
% 
% See also: List related files here

% $Author: suwartad $	$Date: 2017/05/30 16:46:04 $	$Revision: 0.1 $
% Copyright: Process Control Group - NTNU Trondheim 2017

%% setup parameters
clear problem;
sym param;
%isComplete      = 0;
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

etamax = currentEta;

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
success = 0;
testrepeat = 0;
fail    = 0;

% main loop
while t < 1

    numIter = numIter + 1;
    %% compute Eta from previous step
    clear p0 problem;
    p0      = (1 - t)*paramInit + t*paramFinal;
    problem = problemName(p0);
    [objectiveFunctionValue,gradientObjective] = problem.obj(primalInit,p0);
    [constraint,jacobianConstraint]            = problem.cons(primalInit,p0);
    [currentEta,~]                             = calculateEta(problem, dualInit, gradientObjective, constraint, jacobianConstraint);
    %[nextEta,~]                             = calculateEta(problem, dualInit, gradientObjective, constraint, jacobianConstraint);
    %% update parameter
    clear p;
    p  = (1 - t - deltaT)*paramInit + (t + deltaT)*paramFinal;
%     % define nextEta HERE !
%     [objectiveFunctionValue,gradientObjective] = problem.obj(primalInit, p);
%     [constraint,jacobianConstraint]            = problem.cons(primalInit, p);
%     [nextEta,~]                             = calculateEta(problem, dualInit, gradientObjective, constraint, jacobianConstraint);
    
    % change stopping criteria
    %while (nextEta > 1e-4) || ( nextEta > currentEta)
    %while( nextEta > max(currentEta,Parameters.etaMax) )
    if( currentEta < 1 )
     
        %% SOLVE CorrectStep
        [deltaXc,deltaYplus] = solveCorrectStep(problem, variables, constraints, objective, p0);
        
        % check dual variable for inequality constraint if there is any less
        % than zero
        yInequality = variables.dual(inequalityIndex) + deltaYplus(inequalityIndex);
        isNegative  = find(yInequality < 0);
        if (isNegative)
            fprintf('SOLVE NLP AGAIN!!! \n');
            numFail = 1;
            keyboard;
        end
        
        %% SOLVE QPPredict
        [deltaXp,deltaYp, exitQP] = solveQPPredict(problem, variables, constraints, deltaXc, deltaT, deltaP, p0);
        
        % set solution delta_x and delta_y = [delta_xp,delta_yp] + [delta_xc,
        % delta_yc]
        deltaX = deltaXc + deltaXp;
        deltaY = deltaYplus + deltaYp;
        
        %% check condition (5.1)
        clear problem;
        problem = problemName(p);
        [~,gradientObjective]           = problem.obj(primalInit+deltaX, p);
        [constraint,jacobianConstraint] = problem.cons(primalInit+deltaX, p);
        [nextEta,Lag]  = calculateEta(problem, dualInit+deltaY, gradientObjective, constraint, jacobianConstraint);
        
        
        %% testing with new stopping criteria
        %while (nextEta > 1e-4) ||  (nextEta > currentEta^1.2 && nextEta > 1e-2) || ( nextEta > currentEta)
        %while (nextEta > 1e-4) || ( nextEta > currentEta)
        %while( nextEta > max(currentEta,Parameters.etaMax) )
        %while ( (nextEta > 1e-6) ||  (nextEta > currentEta^1.2 && nextEta > 1e-2) )
        
        %if( nextEta <= max(currentEta,Parameters.etaMax) )
        %if (nextEta <= 1e-6) ||  (nextEta <= currentEta^1.2 && nextEta < 1e-2)
        %if (nextEta <= Parameters.etaMax) ||  (nextEta <= currentEta^1.4 && nextEta < max(1e-2,etamax) || (nextEta < etamax) )
        %if nextEta <= max((Parameters.etaMax),etamax) ||  (nextEta <= currentEta^1.4 && nextEta < max(1e-2,etamax)) 
        if nextEta <= max((Parameters.etaMax),etamax)
            
            if numFail > Parameters.maxFailure
                % check failing counter
                fprintf('SOLVE NLP AGAIN!!! \n');
                %numFail = 1;
                keyboard;
            end
            
            % compute estimate of active set A
            % currentEta = nextEta;
            dualInit   = dualInit+deltaY;
            primalInit = primalInit+deltaX;
            % update t
            t          = t + deltaT;
            [activeSet, estimateActive,estimateInactive] = estimateActiveSet(dualInit, constraint, equalityIndex, inequalityIndex, currentEta, Parameters);
            
            % check good eta condition, if satisfied increase delta_t
            if (nextEta <= currentEta^1.6 || nextEta <= Parameters.etaMax) && ~testrepeat
                %if (nextEta < currentEta)
                deltaT = min(1-t,deltaT/Parameters.alpha);
            else
                deltaT = min(1-t,deltaT);
            end
            
            testrepeat = 0;
            fail       = 0;
            currentEta = nextEta;
            
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
            [dualLP,exitLP]      = solveJumpLP(problem, objective, variables, constraints, deltaT, deltaP, Lag);
            
            linobj = problem.dcdp(activeSet,:)*deltaT;
            
            % define A+ set
            %if exitLP > 0
            if exitLP==1 && (linobj'*dualLP < linobj'*dualInit(activeSet)-sqrt(eps)) 
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
                success = 1;
            end
            
            %currentEta = nextEta;
            primal     = [primal;primalInit'];
            dual       = [dual;dualInit'];
            deltaTc    = [deltaTc;deltaT];
            Eta        = [Eta;currentEta];
            tData      = [tData;t];
            if nextEta > etamax
                etamax = nextEta;
            end
        else
            % decrease deltaT
            deltaT  = deltaT * Parameters.alpha;
            if fail == 0
                testrepeat = 1;
            end
            fail = fail+1;
        end
        
        
        
    else
        % decrease deltaT
        deltaT  = deltaT * Parameters.alpha;
    end

%     % check good eta condition, if satisfied increase delta_t
%     if (nextEta < currentEta^(1+Parameters.gamma))
%     %if (nextEta < currentEta)
%         deltaT = min(1-t,deltaT/Parameters.alpha);
%     else
%         deltaT = min(1-t,deltaT);
%     end

    % heuristic to prevent small Eta
    if nextEta < 1e-8
        deltaT = min(1-t,deltaT/Parameters.alpha);
    end
    
    if success == 1
        currentEta = nextEta;
    end


    % print out every iteration ...
    %fprintf('Iter = %6g    Eta = %5.2e    t =  %3.1e    deltaT = %3.1e    QPflag=%d \n', numIter, currentEta, t, deltaT, exitQP);
    fprintf('Iter = %6g    Eta = %5.2e    t =  %3.1e    deltaT = %3.1e \n', numIter, currentEta, t, deltaT);
end

%% PLOT PRIMAL, DUAL, ETA, and deltaT !!!
keyboard;

end

function Parameters = setParameters()
Parameters.gamma        = 0.7;
Parameters.alpha        = 0.6;
%Parameters.alpha        = 0.9;
Parameters.maxIteration = 50000;
Parameters.initDeltaT   = 0.05;
%Parameters.initDeltaT   = 0.1;
Parameters.maxFailure   = 50;
%Parameters.maxFailure   = 40;
Parameters.optValue     = 1e-6;
Parameters.etaMin       = 1e-6;
%Parameters.etaMax       = 1e-2;
%Parameters.etaMax       = 1e-4;
Parameters.etaMax       = 0.05;

end

function [eta,Lagrangian] = calculateEta(problem, dual, grad, cons, Jac)
% compute Lagrangian
Lagrangian = grad - Jac'*dual;

% equality constraints
equalityConstraints   = cons(1:problem.me);
inequalityConstraints = cons(problem.me+1:end);

% % bigger value between inequality constraints and dual variables
% numInequalityConstraints = size(inequalityConstraints,1);
% minInequalityAndDual     = zeros(numInequalityConstraints,1);
% correspondingDual        = dual(problem.me+1:end);
% for i=1:numInequalityConstraints
%     minInequalityAndDual(i,1) = min(inequalityConstraints(i),correspondingDual(i));
% end
% 
% % construct vector of Eta
% vectorEta = [Lagrangian;equalityConstraints;minInequalityAndDual];
% % Eta is infinite norm of that vector
% eta       = norm(vectorEta,inf);

cmp      = max(0,dual(problem.me+1:end).*min(inequalityConstraints,1));
vioL     = min(0,inequalityConstraints);
vioE     = equalityConstraints;


maxComp  = norm( cmp ,Inf);
maxViol  = norm(    [ vioL; vioE ] ,Inf);
eta = .1*max(maxComp,max(maxViol,norm(Lagrangian,Inf)));

end

function [activeSet, estimateActive, estimateInactive]  = estimateActiveSet(dualInit, constraint, equalityIndex, inequalityIndex, Eta, Parameters)
estimateActive   = find(abs(dualInit) > Eta^Parameters.gamma);
estimateActive   = setdiff(estimateActive,equalityIndex);
estimateInactive = setdiff(inequalityIndex ,estimateActive);
activeConstraint = union(equalityIndex,find(constraint < Eta^(Parameters.gamma)));
activeSet        = union(estimateActive,activeConstraint);
end

function [deltaXc,deltaYc] = solveCorrectStep(problem, variables, constraints, objective, p)
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

numAct  = size([J(eqs,:);J(epsA,:)],1);
lhs     = [H                      -[J(eqs,:);J(epsA,:)]'; ...
          [J(eqs,:);J(epsA,:)]   zeros(numAct,numAct)];
rhs     = -[g-J'*y;[c(eqs,:);c(epsA,:)]];
solCorrectStep = lhs\rhs;
deltaXc = solCorrectStep(1:n);
yCS     = solCorrectStep(n+1:end);
deltaYc = zeros(m,1);
deltaYc([eqs;epsA]) = yCS;
end

function [deltaXp,deltaYp, exitflag] = solveQPPredict(problem, variables, constraints, deltaXc, deltaT, deltaP, p)

n       = problem.n;
m       = problem.m;
me      = problem.me;
eqs     = 1:me;
H       = problem.hess(variables.primal,variables.dual,p);
%H       = problem.hess(variables.primal,variables.dual);
HC      = problem.chess(variables.primal);
dcdp    = problem.dcdp;
J       = constraints.jacobian;
% g       = objective.gradient;
% y       = variables.dual;
c       = constraints.value;
epsA    = constraints.eqInd;
epsF    = constraints.inInd;

JacHessianXc = J;
% addition of Jacobian of constraint PLUS Hessian of constraint MULTIPLY BY
% deltaXc
for i=1:m
    JacHessianXc(i,:) = JacHessianXc(i,:) + (squeeze(HC(:,:,i))*deltaXc)';
end
A      = -JacHessianXc(epsF,:);
%b      = dcdp(epsF,:)*(deltaT*deltaP);
b      = c(epsF) + JacHessianXc(epsF,:)*deltaXc + dcdp(epsF,:)*(deltaT*deltaP);
%b      = c(epsF) + dcdp(epsF,:)*(deltaT*deltaP);
Aeq    = [JacHessianXc(eqs,:);JacHessianXc(epsA,:)];
beq    = [-dcdp(eqs,:);-dcdp(epsA,:)]*(deltaT*deltaP);
grad   = zeros(n,1);
option = optimset('Display','off','Algorithm','active-set');
%option = optimset('Display','off');
[deltaXp,~,exitflag,~,lambda] = quadprog(H,grad,A,b,Aeq,beq,[],[], variables.primal, option);
deltaYp             = zeros(m,1);
deltaYp(epsF)       = (lambda.ineqlin(1:(length(epsF))));
deltaYp([eqs;epsA]) = -lambda.eqlin;
% NEED TO CHECK exitflag HERE !!!
% if exitflag < 0
%     deltaXp = zeros(n,1);
%     deltaYp = zeros(m,1);
% else
%     deltaYp             = zeros(m,1);
%     deltaYp(epsF)       = (lambda.ineqlin(1:(length(epsF))));
%     deltaYp([eqs;epsA]) = -lambda.eqlin;
% end
end

function [y,exitflag] = solveJumpLP(problem, objective, variables, constraints, deltaT, deltaP, z)

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
%option   = optimoptions('linprog','Algorithm','dual-simplex','Display','off');
%option   = optimoptions('linprog','Algorithm','simplex','Display','off');
option   = optimoptions('linprog','Algorithm','simplex','Display','iter');
%option   = optimoptions('linprog','Algorithm','active-set','Display','iter');
[y,~,exitflag] = linprog(f, A, b, [],[],lb,ub,variables.dual(Act),option );
%[y,~,exitflag] = linprog(f, A, b, [],[],lb,ub,[],option );
end

