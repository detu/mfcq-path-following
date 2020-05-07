function [Tall, xmeasureAll, uAll, ObjVal, runtime] = spfNmpcCstr(optProblem, system, mpciterations, N, T, tmeasure, xmeasure, u0, scr, varargin)
%SPFNMPCCSTR Summary of this function goes here
% 
% [OUTPUTARGS] = SPFNMPCCSTR(INPUTARGS) Explain usage here
% 
% Examples: 
% 
% Provide sample usage code here
% 
% See also: List related files here

% $Author: suwartad $	$Date: 2018/05/02 23:47:22 $	$Revision: 0.1 $
% Copyright: Process Control Group - NTNU Trondheim 2018

% TO DO:
% 1. Include parameter (price) as measurement data
% 2. Implement as-msNMPC that consists of:
% 2a. Background step (parallel NLP solver; remove read/write file; use
% MATLAB parfor; avoid global variable)
% 2b. Online step: use path-following algorithm

% NEW TO DO:
% RUN ideal-msNMPC at the first iteration in order to hold NAC.
% FROM second iteration run ad-msNMPC.

import casadi.*
%% ITERATE ON MPC RUNNING
% outer loop of the scenario
global nx k pQ pC scenario;
k       = 1.2;
z1      = xmeasure;
t0      = tmeasure;
x0      = xmeasure;
% data initialization
Tall    = [];
Xall    = zeros(mpciterations,size(xmeasure,1));
Uall    = zeros(mpciterations,size(u0,1));
xmeasureAll = [];
uAll        = [];
runtime     = [];
u_pf_opt    = [];
x_pf_opt    = [];

% % Needed for computing Hessian of the constraints
% flagDt  = 1;  

% Initial price information
% pQ = 5;
% pC = 5;
pC = 2;
pQ = 0.5;

% % Start of the NMPC iteration
% global flagDt mpciter;
mpciter = 1;
% flagDt  = 1;  % define 1 one time!
% etaData = [];
% nABC    = [];
while (mpciter <= mpciterations)
    %% printing iteration
    fprintf('-----------------------------\n');
    fprintf('MPC iteration: %d\n', mpciter);
    
    %% FIRST ITERATION CALL ideal-msNMPC
    %if mpciter == 1 || mpciter == 2 || mpciter == 3
    %if mpciter == 1 || mpciter == 2 || mpciter == 25 || mpciter == 26 || mpciter == 80 || mpciter == 81
    if mpciter == 1 || mpciter == 25 || mpciter == 80
        % call sNmpcCstr
        [t0, x0, pQ, pC] = measureInitialValueParam ( tmeasure, xmeasure, mpciter, pQ, pC );
        [u_pf_opt,t0,xmeasure,elapsedqp] = msNmpcCstr(optProblem, N, u0, tmeasure, xmeasure, scr);
    else
        % call advanced-step-multistage-NMPC
        %% BACKGROUND STEP
        %[t0_measure, x0_measure, pQ, pC] = measureInitialValueParam ( tmeasure, xmeasure, mpciter, pQ, pC );
        [t0, x0, pQ, pC] = measureInitialValueParam ( tmeasure, xmeasure, mpciter, pQ, pC );
        
        % ITERATE ON THE SCENARIO TREES (should be parallel, let try serial
        % first)
        %     primalNLP = zeros(452,4);  % hardcode!
        %     dualNLP   = zeros(452,4);
        for i=1:scr.numScr
            % call MPC with v0=uk and z0=xk
            scenario = i;
            [primalNLP(:,i), dualNLP(:,i), lb, ub, ~, params] = solveMsOCP(optProblem, system, N, t0, x0, u0, T, mpciter, u_pf_opt, x_pf_opt, z1);
            %[primalNLP(:,i), dualNLP(:,i), lb, ub, ~, params] = solveMsOCP(optProblem, system, N, t0_measure, x0_measure, u0, T, mpciter, u_pf_opt, x_pf_opt, z1);
        end
        
        % BEFORE SOLVING NLP-SENSITIVITY CHOOSE THE CLOSEST SCENARIO FROM MEASUREMENT
        distance = zeros(scr.numScr,1);
        for i=1:scr.numScr
            %distance(i) = norm(primalNLP(i,8:9) - x0,2); % hardcode position 8 and 9
            %distance(i) = norm(primalNLP(i,3:4) - x0,2);
            %distance(i) = norm(primalNLP(1:2,i) - x0,2);
            %distance(i) = norm(primalNLP(4:5,i) - xmeasure,2);   % without
            %distance(i) = norm(primalNLP(7:8,i) - xmeasure,2);
            %distance(i) = norm(primalNLP(31:32,i) - xmeasure,2);
            distance(i) = norm(primalNLP(22:23,i) - xmeasure,2);
            %distance(i) = norm(primalNLP(7:8,i) - xmeasure,2);
        end
        [~,index] = min(distance);
        fprintf('index number = %d \n', index);
        pfNLP     = primalNLP(:,index);
        dfNLP     = dualNLP(:,index);
        %% ONLINE STEP
        
        %  Obtain new measurement
        %[t0_measure, x0_measure, pQ, pC] = measureInitialValueParam ( tmeasure, xmeasure, mpciter, pQ, pC );
        
        % APPLY PATH-FOLLOWING ALGORITHM HERE!
        % re-arrange NLP solutions
        [~, x_nlp_opt] = arrangeOptResults(pfNLP, lb, ub, N);
        
        %p_init  = pfNLP(1:nx);
        %p_final = x0_measure;
        %p_init  = pfNLP(4:5);
        %p_init  = pfNLP(31:32);
        %p_init  = pfNLP(34:35); % with steady-state variable
        p_init  = pfNLP(22:23);
        %p_init  = pfNLP(7:8);
        p_final = x0;
        %p_final = p_init;
        xstart  = pfNLP;
        ystart  = dfNLP;
        
        % choose number of path-following step
        %delta_t = 0.5;
        delta_t = 1;
        
        lb_init = lb;
        ub_init = ub;
        
        % NLP sensitivity (predictor-corrector)
        [primalPF, ~, elapsedqp] = jpredictor_licq_pure_3(@(p)cstrMultiStage(p), p_init, p_final, xstart, ystart, delta_t, lb_init, ub_init, 0, N);
        %[primalPF, ~, elapsedqpAll, etaRecord, numActiveBoundRecord] = pf_pc_mfcq(@(p)cstrMultiStage(p), p_init, p_final, xstart, ystart, delta_t, lb_init, ub_init, 0, N);
        %elapsedqp = elapsedqpAll.total;
        %etaData = [etaData; etaRecord];
        %nABC    = [nABC; numActiveBoundRecord];
        
        [u_pf_opt, x_pf_opt] = arrangeOptResults(primalPF, lb, ub, N);
        z1 = x_pf_opt(1:nx,5); % 4 = (d+1) + 1 (d=number of collocation point)
    end
    
    %% INJECT TO PLANT
    %x0 = x0_measure; 
    %t0 = t0_measure;
    x0 = xmeasure;
    [tmeasure, xmeasure] = applyControl(system, T, t0, x0, u_pf_opt(1,:));     
    
    %% GET CLOSED-LOOP RESPONSE
    ObjVal(mpciter)      = computeObjFuncCstr(u_pf_opt(1,1),xmeasure); % CSTR only
    
    %% COLLECT CLOSED-LOOP DATA
    Tall                = [Tall t0];
    Xall(mpciter+1,:)   = x0';
    Uall(mpciter+1,:)   = u0(:,1);
    xmeasureAll = [xmeasureAll;xmeasure];
    uAll        = [uAll;u_pf_opt(1,1)];
    runtime     = [runtime;elapsedqp];
    
    %% SHIFT CONTROL INPUT
    u0 = shiftHorizon(u_pf_opt(1,:));
    
    mpciter = mpciter+1;
end
xmeasureAll = reshape(xmeasureAll,nx,mpciterations);

end
