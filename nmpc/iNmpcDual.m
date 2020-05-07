function [Tall, xmeasureAll, uAll, ObjVal, primalNLP, dualNLP, params, runtime] = iNmpcDual(optProblem, system, mpciterations, N, T, tmeasure, xmeasure, u0, xGuess, varargin)
%INMPC Summary of this function goes here
% 
% [OUTPUTARGS] = INMPC(INPUTARGS) Explain usage here
% 
% Ideal NMPC (single step)
%
% Examples: 
% 
% Provide sample usage code here
% 
% See also: List related files here

% $Author: suwartad $	$Date: 2016/04/25 19:24:15 $	$Revision: 0.1 $
% Copyright: Process Control Group - NTNU Trondheim 2016


Tall    = [];
Xall    = zeros(mpciterations,size(xmeasure,1));
Uall    = zeros(mpciterations,size(u0,1));

xmeasureAll = [];
uAll        = [];
runtime     = [];
x_nlp_opt   = [];
u_nlp_opt   = [];
mpciter     = 1;

% storing primal and dual results
primalNLP = [];
dualNLP   = [];


load noise1pct.mat;

global nx;
%xc = 10;
xc = 1;

while(mpciter <= mpciterations)
    
    fprintf('-----------------------------\n');
    fprintf('MPC iteration: %d\n', mpciter);
    
    % obtain new initial value
    [t0, x0] = measureInitialValue ( tmeasure, xmeasure );
    
    % scale bottom concentration
    x0(1) = xc*x0(1);
    
    % mesti cut buat x0 juga!
    x0 = max(min(x0,1.0),0); % restrict to boundaries
    %if x0(1,1) > 0.1; x0(1,1) = 0.1; end
    if x0(1,1) > xc/10; x0(1,1) = xc/10; end        %scaling on bottom concentration
    if x0(84,1) > 0.75; x0(84,1) = 0.75; end
% %      if x0(84,1) < 0.74; x0(84,1) = 0.74; end
% %     if x0(43,1) > 0.75; x0(43,1) = 0.75; end
% %     if x0(43,1) < 0.49; x0(43,1) = 0.49; end
% % %     if x0(83,1) > 0.505; x0(83,1) = 0.505; end
% % %     if x0(83,1) < 0.495; x0(83,1) = 0.495; end
%     if x0(84,1) < 0.7; x0(84,1) = 0.7; end
%     %if x0(84,1) < 0.65; x0(84,1) = 0.65; end        %with measurement noise
%      if x0(84,1) > 0.75; x0(84,1) = 0.75; end
%     if x0(43,1) > 0.6; x0(43,1) = 0.6; end
%     %if x0(43,1) < 0.5; x0(43,1) = 0.5; end
%     if x0(43,1) < 0.4; x0(43,1) = 0.4; end          %with measurement noise
% %     if x0(43,1) > 0.6; x0(43,1) = 0.6; end
% %     if x0(43,1) < 0.5; x0(43,1) = 0.5; end
%      if x0(43,1) > 0.5; x0(43,1) = 0.5; end
%      if x0(43,1) < 0.5; x0(43,1) = 0.5; end
%      if x0(83,1) > 0.5; x0(83,1) = 0.5; end
%      if x0(83,1) < 0.5; x0(83,1) = 0.5; end

    
    % add measurement error to x0
    holdupNoise        = noise(:,mpciter);
    concentrationNoise = zeros(42,1);
    measNoise          = [concentrationNoise;holdupNoise];
    x0_measure         =  x0 + 0*measNoise;    % without noise
    %x0_measure         =  x0 + measNoise;
    
%     % because of level controller at bottom and top, remove the noises
%     measNoise(43)      = 0;
%     measNoise(83)      = 0;
%     x0_measure         =  x0 + measNoise;
    
    
%      dummyNoise         = noise(1:2,mpciter);   % SHOULD GENERATE OWN NOISE !!!
%      x0_measure         =  x0 + dummyNoise;
    
    % check constraints on boundary
    x0_measure = max(min(x0_measure,1.0),0); % restrict to boundaries
% %     %if x0_measure(1,1) > 0.1; x0_measure(1,1) = 0.1; end
     if x0_measure(1,1) > xc/10; x0_measure(1,1) = xc/10; end    %scaling on bottom concentration
     if x0_measure(84,1) > 0.75; x0_measure(84,1) = 0.75; end
%     if x0_measure(84,1) < 0.7; x0_measure(84,1) = 0.7; end
%     %if x0_measure(84,1) < 0.65; x0_measure(84,1) = 0.65; end    %with measurement noise
%     %if x0_measure(84,1) > 0.7; x0_measure(84,1) = 0.7; end
%     %if x0_measure(84,1) > 0.72; x0_measure(84,1) = 0.72; end
%     if x0_measure(43,1) > 0.6; x0_measure(43,1) = 0.6; end
%     %if x0_measure(43,1) < 0.5; x0_measure(43,1) = 0.5; end
%     if x0_measure(43,1) < 0.4; x0_measure(43,1) = 0.4; end     %with measurement noise
% %     if x0_measure(43,1) > 0.6; x0_measure(43,1) = 0.6; end
% %     if x0_measure(43,1) < 0.5; x0_measure(43,1) = 0.5; end
%      if x0_measure(43,1) > 0.5; x0_measure(43,1) = 0.5; end
%      if x0_measure(43,1) < 0.5; x0_measure(43,1) = 0.5; end
%      if x0_measure(83,1) > 0.5; x0_measure(83,1) = 0.5; end
%      if x0_measure(83,1) < 0.5; x0_measure(83,1) = 0.5; end

     

    % ideal NMPC:
    %[primalNLP, ~, lb, ub, ~, params, elapsedtime] = solveOptimalControlProblem(optProblem, system, N, t0, x0, u0, T, mpciter, u_nlp_opt, x_nlp_opt, x0_measure);
    [primalData, dualData, lb, ub, ~, params, elapsedtime] = solveOptimalControlProblem(optProblem, system, N, t0, x0, u0, T, mpciter, u_nlp_opt, x_nlp_opt, x0_measure, xGuess);
    %primalNLP = [primalNLP primalData];
    %dualNLP   = [dualNLP dualData];
    
    %save directMultipleShooting.mat primalNLP dualData;
    
    % re-arrange NLP results
    [u_nlp_opt, x_nlp_opt] = plotStatesN(primalData, lb, ub, N);
    %[u_nlp_opt, x_nlp_opt] = plotStatesSoftConstraint(primalNLP, N);
    %[u_nlp_opt, x_nlp_opt] = plotStatesN(primalData, lb, ub, N);
    %[u_nlp_opt, x_nlp_opt] = plotStatesDirectMS(primalNLP, N);
    primalNLP{mpciter}.u_nlp_opt = u_nlp_opt;
    primalNLP{mpciter}.x_nlp_opt  = x_nlp_opt;
%     xGuess = x_nlp_opt;
%     xGuess(:,1) = [];
    
    % scale bottom concentration
    x_nlp_opt(1,:) = x_nlp_opt(1,:)./xc;
    
    % take the Lagrangian multiplier for bound constraint
    [dual_u, dual_x] = plotStatesN(dualData, lb, ub, N);
    %[dual_u, dual_x] = plotStatesSoftConstraint(dualData, N);
    %[dual_u, dual_x] = plotStatesDirectMS(dualData, N);
    dualNLP{mpciter}.dualU = dual_u;
    dualNLP{mpciter}.dualX = dual_x;
    
%     if mod(mpciter,10) == 0
%         keyboard;
%     end
    
    % save open loop solution for error computation
    %iNmpcData(mpciter).z = x_nlp_opt;

    z1 = x_nlp_opt(1:nx,5);  % 5 = (d+1) + 1 (d=number of collocation point)
    
    
    % record information
    Tall                = [Tall t0];
    Xall(mpciter+1,:)   = x0';
    Uall(mpciter+1,:)   = u0(:,1);
    
    
    % Apply control to process with optimized control from path-following
    % algorithm. 
    %x0 = xmeasure;  % from the online step 
    x0 = x0_measure; % NEW CHANGE 28.09.2017
    x0(1) = x0(1) / xc; % scale bottom concentration
    [tmeasure, xmeasure] = applyControl(system, T, t0, x0, u_nlp_opt);
    
    ObjVal(mpciter) = computeObjectiveFunctionValues(u_nlp_opt(:,1),xmeasure); % 
    %ObjVal(mpciter) = mpciter; % DUMMY!
    %ObjVal(mpciter) = computeObjFuncCstr(u_nlp_opt(:,1),xmeasure); % CSTR only
    
    % store output variables
    xmeasureAll = [xmeasureAll;xmeasure];
    uAll        = [uAll;u_nlp_opt(:,1)];   
    runtime     = [runtime;elapsedtime];
    
    % prepare restart
    u0 = shiftHorizon(u_nlp_opt);
    
    mpciter = mpciter+1;
end
xmeasureAll = reshape(xmeasureAll,nx,mpciterations);    
%save iNmpcData.mat iNmpcData;
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

function [u, lamda, lbw, ubw, objVal, params, elapsednlp] = solveOptimalControlProblem(optProblem, system, N, t0, x0, u0, T, mpciter, u_nlp, x_nlp, x0_measure, xGuess)

%     load Xinit31.mat;
%     x0 = Xinit31(1:84);
%     
     import casadi.*
%     % call ODE15s N-times initial guess in optimization
%     x(1,:) = x0';
%     for k=1:N
%         x(k+1,:) = x0';    % ONE-TIME SIMULATION !
%     end

    x = xGuess;
%     load Xopt31.mat;
%     %load Xopt3175.mat;
%     x = Xopt;

%     x = repmat(x0,120,1);
%     x = reshape(x,84,120);

%     load Xopt3175Horizon90.mat
%     load Xopt3164Horizon90.mat;
%     x = Xopt;

    [J,g,w0,w,lbg,ubg,lbw,ubw,params] = optProblem(x, u0, N, x0_measure);
    prob = struct('f', J, 'x', vertcat(w{:}), 'g', vertcat(g{:}));
    options = struct;
    %options.ipopt.tol                = 1e-12;
    %options.ipopt.constr_viol_tol    = 1e-6;
    %options.ipopt.constr_viol_tol    = 1e-12;  
    %options.ipopt.dual_inf_tol       = 1e-12;
%     options.ipopt.mu_strategy        = 'adaptive';
%     options.ipopt.compl_inf_tol      = 1e-12;
%     options.ipopt.acceptable_tol     = 1e-12;  
    %options.ipopt.nlp_scaling_method  = 'gradient-based'; 
    %options.ipopt.nlp_scaling_min_value = 1e-5;
%     options.ipopt.least_square_init_duals    = 'yes';
%     options.ipopt.adaptive_mu_globalization  = 'kkt-error';
    
    solver = nlpsol('solver', 'ipopt', prob, options);
    %solver = nlpsol('solver', 'scpgen', prob, options);
    %solver = nlpsol('solver', 'sqpmethod', prob, options);
    
    %solver = nlpsol('solver', 'knitro', prob, options);

    
    % Solve the NLP
    startnlp   = tic;
    sol        = solver('x0', w0, 'lbx', lbw, 'ubx', ubw, 'lbg', lbg, 'ubg', ubg);
    elapsednlp = toc(startnlp);
    %elapsednlp = solver.stats.t_wall_mainloop;
    fprintf('IPOPT solver runtime = %f\n',elapsednlp);
    success = strcmp(solver.stats.return_status,'Infeasible_Problem_Detected');
    if (success)
        keyboard;
    end

    u      = full(sol.x);
    lamdaG  = full(sol.lam_g);
    lamda  = full(sol.lam_x);
    objVal = full(sol.f);
end

function [x, t_intermediate, x_intermediate] = dynamic(system, T, t0, x0, u0)
    x = system(t0, x0, u0, T);
    x_intermediate = [x0; x];
    t_intermediate = [t0, t0+T];
end

