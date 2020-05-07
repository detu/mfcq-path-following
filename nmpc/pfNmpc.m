function [Tall, xmeasureAll, uAll, ObjVal, primalPF, params, runtime, etaData, nABC, activeChange] = pfNmpc(optProblem, system, mpciterations, N, T, tmeasure, xmeasure, u0, varargin)
%PFNMPC Summary of this function goes here
% 
% Path-following based Nonlinear Model Predictive Control
%
% [OUTPUTARGS] = PFNMPC(INPUTARGS) Explain usage here
% 
% Examples: 
% 
% Provide sample usage code here
% 
% See also: List related files here

% $Author: suwartad $	$Date: 2016/04/06 20:39:51 $	$Revision: 0.1 $
% Copyright: Process Control Group - NTNU Trondheim 2016

% dimension of state x and control u
nx     = numel(xmeasure);
nu     = size(u0,1);
Tall   = [];
Xall   = zeros(mpciterations,size(xmeasure,1));
Uall   = zeros(mpciterations,size(u0,1));

xmeasureAll = [];
uAll        = [];
runtime     = [];

u_pf_opt    = [];
x_pf_opt    = [];

% measuring active set during the course of pf-NMPC iteration
activeChange = zeros(mpciterations,1);
countActive  = 0;

load noise1pct.mat;
z1 = xmeasure;
% Start of the NMPC iteration
global flagDt mpciter;
mpciter = 1;
flagDt  = 1;  % define 1 one time!

etaData = [];
nABC    = [];
while(mpciter <= mpciterations)
    
    fprintf('-----------------------------\n');
    fprintf('MPC iteration: %d\n', mpciter);

    %   Obtain new initial value
    [t0, x0] = measureInitialValue ( tmeasure, xmeasure );
    
%     % mesti cut buat x0 juga!
    x0 = max(min(x0,1.0),0); % restrict to boundaries
    if x0(1,1) > 0.1; x0(1,1) = 0.1; end
    if x0(84,1) > 0.75; x0(84,1) = 0.75; end
    if x0(84,1) < 0.74; x0(84,1) = 0.74; end
    if x0(43,1) > 0.5; x0(43,1) = 0.5; end
    if x0(43,1) < 0.5; x0(43,1) = 0.5; end
    
%     if x0(2,1) < 0.49; x0(2,1) = 0.49; end 
    
    % add measurement noise
    holdupNoise        = noise(:,mpciter);
    concentrationNoise = zeros(42,1);
    measNoise          = [concentrationNoise;holdupNoise];
    x0_measure         =  x0 + 0*measNoise;  % without noise
    %x0_measure         =  x0 + measNoise;
    
    % check constraints on boundary
    x0_measure = max(min(x0_measure,1.0),0); % restrict to boundaries
    if x0_measure(1,1) > 0.1; x0_measure(1,1) = 0.1; end
    if x0_measure(84,1) > 0.75; x0_measure(84,1) = 0.75; end
    if x0_measure(84,1) < 0.74; x0_measure(84,1) = 0.74; end
    if x0_measure(43,1) > 0.5; x0_measure(43,1) = 0.5; end
    if x0_measure(43,1) < 0.5; x0_measure(43,1) = 0.5; end
    
    
    z1 = max(min(z1,1.0),0); % restrict to boundaries
    if z1(1,1) > 0.1; z1(1,1) = 0.1; end
    if z1(84,1) > 0.75; z1(84,1) = 0.75; end
    if z1(84,1) < 0.74; z1(84,1) = 0.74; end
    if z1(43,1) > 0.5; z1(43,1) = 0.5; end
    if z1(43,1) < 0.5; z1(43,1) = 0.5; end



%      dummyNoise         = noise(1:2,mpciter);   % SHOULD GENERATE OWN NOISE !!!
%      x0_measure         =  x0 + dummyNoise;
%      if x0_measure(2,1) < 0.49; x0_measure(2,1) = 0.49; end


    % advanced-step NMPC:
    [primalNLP, dualNLP, lb, ub, ~, params] = solveOptimalControlProblem(optProblem, system, N, t0, x0, u0, T, mpciter, u_pf_opt, x_pf_opt, z1);
    
    % re-arrange NLP solutions
    [~, x_nlp_opt] = plotStatesN(primalNLP, lb, ub, N);
    
    p_init  = primalNLP(1:nx);
    %p_init  = x0;
    p_final = x0_measure; 
    xstart  = primalNLP;
    ystart  = dualNLP;
    
    % check if there is any active-set change here !
    if (~isequal(p_init(1),0.1)) || (~isequal(p_final(1),0.1)) 
        if (p_init(1) < 0.1) && (p_final(1) < 0.1)
            countActive = 0;
        else
            countActive = countActive + 1;
        end
    end
    if (~isequal(p_init(84),0.7)) || (~isequal(p_final(84),0.7))
        if (p_init(84) < 0.7) && (p_final(84) < 0.7)
            countActive = countActive + 0;
        else
            countActive = countActive + 1;
        end
    end
    activeChange(mpciter,1) = countActive;
    countActive             = 0;
    
    % choose number of path-following step
    %delta_t = 0.5;   % initial deltaT OK with eta_max = 0.88 
    %delta_t = 0.2;
    %delta_t = 0.25;
    %delta_t = 0.1;
    %delta_t = 0.05;
    delta_t = 1;
    
    lb_init = lb;
    ub_init = ub;
    
    % NLP sensitivity (predictor-corrector)
    % with earlier CasADi 3.1.0
    %[primalPF, ~, elapsedqp, etaRecord, numActiveBoundRecord] = pf_pc_mfcq(@(p)distColACstr_mfcq(p), p_init, p_final, xstart, ystart, delta_t, lb_init, ub_init, 0, N);
    % with CasADi version 3.4.5
    [primalPF, ~, elapsedqp, etaRecord, numActiveBoundRecord] = pf_pc_mfcq(@(p)distColACstr_mfcqNew(p), p_init, p_final, xstart, ystart, delta_t, lb_init, ub_init, 0, N);
    
    % Implement new path-following without corrector step:
    % 1. Predictor-corector QP
    % 2. LP
    %[primalPF, ~, elapsedqp, etaRecord, numActiveBoundRecord] = pf_pc_qp_mfcq(@(p)distColACstr_mfcqNew(p), p_init, p_final, xstart, ystart, delta_t, lb_init, ub_init, 0, N); 

    etaData = [etaData; etaRecord];
    nABC    = [nABC; numActiveBoundRecord];
     
    %[u_nlp_opt, x_nlp_opt] = plotStatesN(primalNLP, lb, ub, N);  % distillation column plot
    [u_pf_opt, x_pf_opt] = plotStatesN(primalPF, lb, ub, N);
    z1 = x_pf_opt(1:nx,5); % 4 = (d+1) + 1 (d=number of collocation point)

    % store output variables
    Tall   = [Tall t0];
    Xall(mpciter+1,:)   = x0';
    Uall(mpciter+1,:)   = u0(:,1);
    
    % Apply control to process with optimized control from path-following
    % algorithm. 
    %x0                   = xmeasure;  % from the online step 
    x0                   = x0_measure; % NEW CHANGE 28.09.2017
    [tmeasure, xmeasure] = applyControl(system, T, t0, x0, u_pf_opt);
    
    ObjVal(mpciter) = computeObjectiveFunctionValues(u_pf_opt(:,1),xmeasure); % USING ACTUAL STATE!
    %ObjVal(mpciter) = computeObjFuncCstr(u_pf_opt(:,1),xmeasure); % CSTR only
    
    
    % collect variables
    xmeasureAll = [xmeasureAll;xmeasure];
    uAll        = [uAll;u_pf_opt(:,1)];   
    runtime     = [runtime;elapsedqp];
    
    % prepare restart: u_new is from path-following output !
    u0      = shiftHorizon(u_pf_opt);
    mpciter = mpciter+1;
end

xmeasureAll = reshape(xmeasureAll,nx,mpciterations);  

end

function [t0, x0] = measureInitialValue ( tmeasure, xmeasure )
    t0 = tmeasure;
    x0 = xmeasure;
end

function [tapplied, xapplied] = applyControl(system, T, t0, x0, u0)
    xapplied = dynamic(system, T, t0, x0, u0(:,1));
    %xapplied = dynamic(system, T, t0, x0, round(u0(:,1),1));
    tapplied = t0+T;
end

function u0 = shiftHorizon(u)
    u0 = [u(:,2:size(u,2)) u(:,size(u,2))]; 
end

function [u, lamda, lbw, ubw, objVal, params] = solveOptimalControlProblem(optProblem, system, N, t0, x0, u0, T, mpciter, u_nlp, x_nlp, x0_measure)

    import casadi.* 
    global ns;
    
    % call ODE15s N-times initial guess in optimization
    x(1,:) = x0';
    for k=1:N
        x(k+1,:) = x0';    % ONE-TIME SIMULATION !
    end
    
    [J,g,w0,w,lbg,ubg,lbw,ubw,params] = optProblem(x, u0, N, x0_measure);
   
    prob    = struct('f', J, 'x', vertcat(w{:}), 'g', vertcat(g{:}));
    options = struct;
    options.ipopt.tol                = 1e-12;
    options.ipopt.constr_viol_tol    = 1e-10; 
    solver = nlpsol('solver', 'ipopt', prob, options);

    % Solve the NLP
    startnlp = tic;
    sol   = solver('x0', w0, 'lbx', lbw, 'ubx', ubw, 'lbg', lbg, 'ubg', ubg);
    elapsednlp = toc(startnlp);
    fprintf('IPOPT solver runtime = %f\n',elapsednlp);

    u           = full(sol.x);
    lamda.lam_g = full(sol.lam_g);
    lamda.lam_x = full(sol.lam_x);
    objVal      = full(sol.f);
end

function [x, t_intermediate, x_intermediate] = dynamic(system, T, t0, x0, u0)
    x = system(t0, x0, u0, T);
    x_intermediate = [x0; x];
    t_intermediate = [t0, t0+T];
end
