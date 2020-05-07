function [Tall, xmeasureAll, uAll, ObjVal, runtime] = sNmpcCstr(optProblem, system, mpciterations, N, T, tmeasure, xmeasure, u0, scr, varargin)
%SNMPCCSTR Summary of this function goes here
%
% Scenario tree NMPC
%
% [OUTPUTARGS] = SNMPCCSTR(INPUTARGS) Explain usage here
%
% Examples:
%
% Provide sample usage code here
%
% See also: List related files here

% $Author: suwartad $	$Date: 2018/04/18 05:31:23 $	$Revision: 0.1 $
% Copyright: Process Control Group - NTNU Trondheim 2018
import casadi.*
%% ITERATE ON MPC RUNNING
% outer loop of the scenario
global nx k xf uf check pQ pC scenario;
xf      = xmeasure;
uf      = u0(1);
mpciter = 1;

% data initialization
Tall    = [];
Xall    = zeros(mpciterations,size(xmeasure,1));
Uall    = zeros(mpciterations,size(u0,1));
xmeasureAll = [];
uAll        = [];
runtime     = [];
x_nlp_opt   = [];
u_nlp_opt   = [];

while (mpciter <= mpciterations)
    %% printing iteration
    fprintf('-----------------------------\n');
    fprintf('MPC iteration: %d\n', mpciter);
    
    %% parameter changes
    % change parameter value
    if mpciter == 1
        check = 0;
        k     = 1.2;
        pC    = 2;
        pQ    = 0.5;
    end
    
    if mpciter == 25
        check = 0;
        pC    = 4;
        pQ    = 0.25;
    end
    
    if mpciter == 80
        check = 0;
        pC    = 1;
        pQ    = 1;
    end
    
    %% ITERATE ON THE SCENARIO TREES
    % prepare optimization variables
    J   = 0;
    g   = {};
    w0  = [];
    w   = {};
    lbg = [];
    ubg = [];
    lbw = [];
    ubw = [];
    nac = [];
    for i=1:scr.numScr
    %for i=1:4
        % call NMPC for each scenario
        % divide into 2 realizations: +- 10%
        
        %% GATHER ALL VARIABLES FROM EACH SCENARIO
        % sNmpcCstr(optProblem, system, mpciterations, N, T, tmeasure, xmeasure, u0, scr, varargin)
        %[Tall, xmeasureAll, uAll, ObjVal, primalNLP, params, runtime] = CollectVariablesForEachScenario(optProblem, system, mpciterations, N, T, u0, tmeasure, xmeasure, i);
        [Ji,gi,w0i,wi,lbgi,ubgi,lbwi,ubwi,t0,x0_measure] = CollectVariablesForEachScenario(optProblem, N, u0, tmeasure, xmeasure, i);
        
        %% collection optimization variables for all realization
        J   = J + Ji;
        %g   = [g;gi];
        g   = {g{:},gi{:}};
        w0  = [w0;w0i];
        w   = {w{:},wi{:}};
        lbg = [lbg;lbgi];
        ubg = [ubg;ubgi];
        lbw = [lbw;lbwi];
        ubw = [ubw;ubwi];
        nac = [nac;wi{3}];  %control input
        %nac = [nac;wi{2}];   % remember control input index!
    end
    
    %% BUILD NON-ANTICIPATIVITY CONSTRAINTS (NAC)
    % ingat ada steady-state constraint JANGAN DIPAKE sebagai NAC!
    % let try robust horizon = 1
    numNac = size(nac,1);
    for i=1:(numNac-1)
        g   = {g{:}, nac(i+1,1)-nac(i,1)};
        lbg = [lbg; 0];
        ubg = [ubg; 0];
    end
    
    
    %% SOLVE NLP FOR ALL REALIZATIONS
    prob = struct('f', J, 'x', vertcat(w{:}), 'g', vertcat(g{:}));
    options = struct;
    options.ipopt.tol                = 1e-12;
    options.ipopt.constr_viol_tol    = 1e-12;
    solver = nlpsol('solver', 'ipopt', prob, options);
    
    % Solve the NLP
    startnlp   = tic;
    sol        = solver('x0', w0, 'lbx', lbw, 'ubx', ubw, 'lbg', lbg, 'ubg', ubg);
    elapsednlp = toc(startnlp);
    fprintf('IPOPT solver runtime = %f\n',elapsednlp);
    success = strcmp(solver.stats.return_status,'Infeasible_Problem_Detected');
    if (success)
        keyboard;
    end
    
    u      = full(sol.x);
    lamda  = full(sol.lam_g);
    objVal = full(sol.f);
    
    %% INJECT TO PLANT
    % reshape optimized variables
    numCol = size(u,1)/scr.numScr;
    u      = reshape(u,numCol,scr.numScr);
    % loop for all scenarios for plotting purposes
    u_nlp_opt  = zeros(scr.numScr,50); % 150 is hardcore
    x_nlp_opt_a = zeros(scr.numScr,201); % idem
    x_nlp_opt_b = zeros(scr.numScr,201);
    for i= 1:scr.numScr
        uopt   = u(:,i);
        [u_nlp_opt(i,:), x_nlp_opt, ~] = arrangeOptResults(uopt, lbw(1:numCol,1), ubw(1:numCol,1), N);
        x_nlp_opt_a(i,:) = x_nlp_opt(1,:);
        x_nlp_opt_b(i,:) = x_nlp_opt(2,:);
    end
    x0 = x0_measure; 
    [tmeasure, xmeasure] = applyControl(system, T, t0, x0, u_nlp_opt(1,:));     
    
    %% GET CLOSED-LOOP RESPONSE
    ObjVal(mpciter)      = computeObjFuncCstr(u_nlp_opt(1,1),xmeasure); % CSTR only
    
    %% COLLECT CLOSED-LOOP DATA
    Tall                = [Tall t0];
    Xall(mpciter+1,:)   = x0';
    Uall(mpciter+1,:)   = u0(:,1);
    xmeasureAll = [xmeasureAll;xmeasure];
    uAll        = [uAll;u_nlp_opt(1,1)];
    runtime     = [runtime;elapsednlp];
    
    %% SHIFT CONTROL INPUT
    u0 = shiftHorizon(u_nlp_opt(1,:));
    
    mpciter = mpciter+1;
end
xmeasureAll = reshape(xmeasureAll,nx,mpciterations);
end

