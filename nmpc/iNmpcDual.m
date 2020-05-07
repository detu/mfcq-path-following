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

    
    % add measurement error to x0
    holdupNoise        = noise(:,mpciter);
    concentrationNoise = zeros(42,1);
    measNoise          = [concentrationNoise;holdupNoise];
    x0_measure         =  x0 + 0*measNoise;    % without noise
    %x0_measure         =  x0 + measNoise;
    
    
    % check constraints on boundary
    x0_measure = max(min(x0_measure,1.0),0); % restrict to boundaries
    %if x0_measure(1,1) > xc/10; x0_measure(1,1) = xc/10; end    %scaling on bottom concentration
    %if x0_measure(84,1) > 0.75; x0_measure(84,1) = 0.75; end
     
    % ideal NMPC:
    [primalData, dualData, lb, ub, ~, params, elapsedtime] = solveOptimalControlProblem(optProblem, system, N, t0, x0, u0, T, mpciter, u_nlp_opt, x_nlp_opt, x0_measure, xGuess);

    
    % re-arrange NLP results
    [u_nlp_opt, x_nlp_opt] = plotStatesN(primalData, lb, ub, N);
    primalNLP{mpciter}.u_nlp_opt = u_nlp_opt;
    primalNLP{mpciter}.x_nlp_opt  = x_nlp_opt;
    
    % scale bottom concentration
    x_nlp_opt(1,:) = x_nlp_opt(1,:)./xc;
    
    % take the Lagrangian multiplier for bound constraint
    [dual_u, dual_x] = plotStatesN(dualData, lb, ub, N);
    dualNLP{mpciter}.dualU = dual_u;
    dualNLP{mpciter}.dualX = dual_x;
   

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
    
    % store output variables
    xmeasureAll = [xmeasureAll;xmeasure];
    uAll        = [uAll;u_nlp_opt(:,1)];   
    runtime     = [runtime;elapsedtime];
    
    % prepare restart
    u0 = shiftHorizon(u_nlp_opt);
    
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
    tapplied = t0+T;
end

function u0 = shiftHorizon(u)
    u0 = [u(:,2:size(u,2)) u(:,size(u,2))];
end

function [u, lamda, lbw, ubw, objVal, params, elapsednlp] = solveOptimalControlProblem(optProblem, system, N, t0, x0, u0, T, mpciter, u_nlp, x_nlp, x0_measure, xGuess)
    
     import casadi.*

    x = xGuess;

    [J,g,w0,w,lbg,ubg,lbw,ubw,params] = optProblem(x, u0, N, x0_measure);
    prob = struct('f', J, 'x', vertcat(w{:}), 'g', vertcat(g{:}));
    options = struct;
    
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
    lamdaG  = full(sol.lam_g);
    lamda  = full(sol.lam_x);
    objVal = full(sol.f);
end

function [x, t_intermediate, x_intermediate] = dynamic(system, T, t0, x0, u0)
    x = system(t0, x0, u0, T);
    x_intermediate = [x0; x];
    t_intermediate = [t0, t0+T];
end

