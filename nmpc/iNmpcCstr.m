function [Tall, xmeasureAll, uAll, ObjVal, primalNLP, params, runtime] = iNmpcCstr(optProblem, system, mpciterations, N, T, tmeasure, xmeasure, u0, varargin)
%INMPCCSTR Summary of this function goes here
% 
% [OUTPUTARGS] = INMPCCSTR(INPUTARGS) Explain usage here
% 
% Examples: 
% 
% Provide sample usage code here
% 
% See also: List related files here

% $Author: suwartad $	$Date: 2018/01/16 22:50:27 $	$Revision: 0.1 $
% Copyright: Process Control Group - NTNU Trondheim 2018

Tall    = [];
Xall    = zeros(mpciterations,size(xmeasure,1));
Uall    = zeros(mpciterations,size(u0,1));

xmeasureAll = [];
uAll        = [];
runtime     = [];
x_nlp_opt   = [];
u_nlp_opt   = [];
mpciter = 1;

load noise1pct.mat;

global nx k xf uf check pQ pC;
xf = xmeasure;
uf = u0(1);

while(mpciter <= mpciterations)
    
    fprintf('-----------------------------\n');
    fprintf('MPC iteration: %d\n', mpciter);
    
    % change parameter value
    if mpciter == 1
        check = 0;
        k  = 1.2;
        pC = 2;
        pQ = 0.5;
    end
    
    if mpciter == 25
        check = 0;
        %k = 0.9;
        pC = 4;
        pQ = 0.25;
    end
    
    if mpciter == 80
        check = 0;
        %k = 2.5;
        pC = 1;
        pQ = 1;
    end
    
    % obtain new initial value
    [t0, x0] = measureInitialValue ( tmeasure, xmeasure );
    
%     % mesti cut buat x0 juga!
    x0 = max(min(x0,1.0),0); % restrict to boundaries
%     if x0(1,1) > 0.1; x0(1,1) = 0.1; end
%     if x0(84,1) > 0.7; x0(84,1) = 0.7; end
%     if x0(2,1) < 0.49; x0(2,1) = 0.49; end
    
    % add measurement error to x0
%     holdupNoise        = noise(:,mpciter);
%     concentrationNoise = zeros(42,1);
%     measNoise          = [concentrationNoise;holdupNoise];
    x0_measure         =  x0;    % without noise
    %x0_measure         =  x0 + measNoise;
%      dummyNoise         = noise(1:2,mpciter);   % SHOULD GENERATE OWN NOISE !!!
%      x0_measure         =  x0 + dummyNoise;
    
%     % check constraints on boundary
    x0_measure = max(min(x0_measure,1.0),0); % restrict to boundaries
%     if x0_measure(1,1) > 0.1; x0_measure(1,1) = 0.1; end
%     if x0_measure(84,1) > 0.7; x0_measure(84,1) = 0.7; end
%     if x0_measure(2,1) < 0.49; x0_measure(2,1) = 0.49; end
     

    % ideal NMPC:
    [primalNLP, ~, lb, ub, ~, params, elapsedtime] = solveOptimalControlProblem(optProblem, system, N, t0, x0, u0, T, mpciter, u_nlp_opt, x_nlp_opt, x0_measure);
    
    % re-arrange NLP results
    %[u_nlp_opt, x_nlp_opt] = plotStatesN(primalNLP, lb, ub, N);  
    [u_nlp_opt, x_nlp_opt, xuSS] = arrangeOptResults(primalNLP, lb, ub, N);
    xuSSt = [x_nlp_opt(end-1);x_nlp_opt(end);u_nlp_opt(end)];
    
    if mpciter == 1 || mpciter == 25 || mpciter == 80
        check = 1;
        xf = [xuSS(1); xuSS(2)];
        uf = xuSS(3);
    end
    
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
    [tmeasure, xmeasure] = applyControl(system, T, t0, x0, u_nlp_opt);
    
    %ObjVal(mpciter) = computeObjectiveFunctionValues(u_nlp_opt(:,1),xmeasure); % 
    %ObjVal(mpciter) = mpciter; % DUMMY!
    ObjVal(mpciter) = computeObjFuncCstr(u_nlp_opt(:,1),xmeasure); % CSTR only
    
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

function [u, lamda, lbw, ubw, objVal, params, elapsednlp] = solveOptimalControlProblem(optProblem, system, N, t0, x0, u0, T, mpciter, u_nlp, x_nlp, x0_measure)

    import casadi.*
    % call ODE15s N-times initial guess in optimization
    x(1,:) = x0';
    for k=1:N
        x(k+1,:) = x0';    % ONE-TIME SIMULATION !
    end

    [J,g,w0,w,lbg,ubg,lbw,ubw,params] = optProblem(x, u0, N, x0_measure);
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
end

function [x, t_intermediate, x_intermediate] = dynamic(system, T, t0, x0, u0)
    x = system(t0, x0, u0, T);
    x_intermediate = [x0; x];
    t_intermediate = [t0, t0+T];
end
