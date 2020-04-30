function NMPCDistCstr
%NMPCDISTCSTR Summary of this function goes here
% 
% NMPC for CSTR + Distillation Column A 
%
% [OUTPUTARGS] = NMPCDISTCSTR(INPUTARGS) Explain usage here
% 
% Examples: 
% 
% Provide sample usage code here
% 
% See also: List related files here

% $Author: suwartad $	$Date: 2016/06/30 11:42:49 $	$Revision: 0.1 $
% Copyright: Process Control Group - NTNU Trondheim 2016

import casadi.* 
format long;

global N;
% number of mpc iteration
mpciterations = 150;
% number of prediction horizon
N             = 45;  
% sampling time
T             = 1;  % [minute]
% initial controls (different initial conditions)
load Xinit30.mat
u0            = Xinit30(85:89);
u0            = repmat(u0,1,N);
% get initial measurement (states) at time T = 0.
tmeasure      = 0.0;
xmeasure      = Xinit30(1:84);

% either call iNMPC 
[~, xmeasureAll, uAll, obj, optRes, params, runtime] = iNmpc(@optProblem, @system, mpciterations, N, T, tmeasure, xmeasure, u0);

% or pf-NMPC
%[~, xmeasureAll_pf, uAll_pf, obj_pf, optRes_pf, params_pf, runtime_pf, etaRecord, numActiveBoundRecord] = pfNmpc(@optProblem, @system, mpciterations, N, T, tmeasure, xmeasure, u0);

keyboard;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Definition of the NMPC functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = system(t, x, u, T)
    
    global uc;
    uc = u;
    [~,x_out] = ode15s('cola_lv_cstr',[t t+T], x);
    lengthx   = size(x_out); 
    y         = x_out(lengthx(1),:)'; 
    
end

function [J,g,w0,w,lbg,ubg,lbw,ubw,params] = optProblem(x, u, N, x0_measure)   %add prediction horizon 
    import casadi.*
    
    % the model
    NT = 41;
    Uf = 0.3;           % Feeding rate F_0
    
    % invoke the model
    [~,state,xdot,inputs] = DistColACstr(Uf);
    f = Function('f',{state,inputs}, {xdot});
    
    % bound constraints
    VB_max = 4.008;
    xB_max = 0.1;
    
    
    % State bounds and initial guess
    x_min =  zeros(84,1);  % try without epsilon here, later put epsilon
    x_max =  ones(84,1);
    
    x_max(1)  = xB_max;
    x_max(84) = 0.75;
    
    % Control bounds
    u_min = [0.1; 0.1; 0.1; 0.1; 0.1];
    u_max = [10; VB_max; 10; 1.0; 1.0];
    
    % compact bound variable for ease of function invocation 
    params.bound.x_min = x_min;
    params.bound.x_max = x_max;
    params.bound.u_min = u_min;
    params.bound.u_max = u_max;
    
    % Construct objective function
    load CstrDistXinit.mat;
    xf    = Xinit(1:84);
    u_opt = Xinit(85:89);
    
    % prices
    pf = 1; 
    pV = 0.02;
    pB = 2; 
    pD = 0;

    
    % compact price variable
    params.price.pf = pf;
    params.price.pV = pV;
    params.price.pB = pB;
    params.price.pD = pD;
    params.price.F_0= Uf;

    % controller gains
    KcB = 10;  
    KcD = 10;
    % Nominal holdups - these are rather small 
    MDs = 0.5; 
    MBs = 0.5;         
    % Nominal flows
    Ds  = 0.5; 
    Bs  = 0.5;
    
    % compact controller gain variable
    params.gain.KcB = KcB;
    params.gain.KcD = KcD;
    params.gain.MDs = MDs;
    params.gain.MBs = MBs;
    params.gain.Ds  = Ds;
    params.gain.Bs  = Bs;
    
    % dimensions
    global nx nu nk d tf ns;
    nx = 84;   % CSTR + Distillation Column A
    nu = 5;    % LT, VB, F, D, B
    nk = 1;
    tf = 1;   % in [minutes]
    h  = tf/nk;
    ns = 0;
    
    % compact model variable
    params.model.NT = NT;
    params.model.f  = f;
    params.model.xdot_val_rf_ss = xf;
    params.model.x  = x;
    params.model.u_opt = u_opt;
    % duplicate u0
    params.model.u  = repmat(u,1,nk); 


    % preparing collocation matrices
    [~,C,D,d] = collocationSetup();
    
    % compact collocation variable
    params.colloc.C = C;
    params.colloc.D = D;
    params.colloc.h = h;

    % start with an empty NLP
    w   = {};      % decision variables contain both control and state variables
    w0  = [];      % initial guess
    lbw = [];      % lower bound for decision variable
    ubw = [];      % upper bound
    J   = 0;       % objective function
    g   = {};      % nonlinear constraint
    lbg = [];      % lower bound for nonlinear constraint
    ubg = [];      % upper bound

    %delta_time = 60; % [minute] convert second to minute
    delta_time = 1;
    alpha = 1;
    beta  = 1;
    gamma = 1;
    
    % compact weight variable
    params.weight.delta_time = delta_time;
    params.weight.alpha      = alpha;
    params.weight.beta       = beta;
    params.weight.gamma      = gamma;
    
    % "Lift" initial conditions
    X0  = MX.sym('X0', nx);
    w   = {w{:}, X0};
    lbw = [lbw; x_min];
    ubw = [ubw; x_max];
    w0  = [w0; x(1,1:nx)'];
    g   = {g{:}, X0 - x0_measure};  % USE MEASUREMENT HERE !
    lbg = [lbg; zeros(nx,1)];
    ubg = [ubg; zeros(nx,1)];

    % formulate the NLP
    Xk = X0;

    
    load Qmax.mat;
    params.Qmax = Qmax;
    
    % HERE SHOULD BE LOOP N-TIMES ACCORDING TO THE NUMBER OF PREDICTION HORIZON
    count  = 2; % counter for state variable as initial guess
    ssoftc = 0;
    for i=1:N
        [J,g,w0,w,lbg,ubg,lbw,ubw,Xk,params,count,ssoftc] = iterateOnPredictionHorizon(Xk, w, w0, lbw, ubw, lbg, ubg, g, J, params, i, count,ssoftc);
    end

    
end


function [J,g,w0,w,lbg,ubg,lbw,ubw,Xk,params,count,ssoftc] = iterateOnPredictionHorizon(Xk, w, w0, lbw, ubw, lbg, ubg, g, J, params, iter, count, ssoftc)

   import casadi.*
   global N;
   % extract compact variables
   x_min = params.bound.x_min;
   x_max = params.bound.x_max;
   u_min = params.bound.u_min;
   u_max = params.bound.u_max;
   
   NT = params.model.NT;
   f  = params.model.f;
   xdot_val_rf_ss = params.model.xdot_val_rf_ss;
   x = params.model.x;
   u = params.model.u;
   u_opt = params.model.u_opt;
   
   pf = params.price.pf;
   pV = params.price.pV;
   pB = params.price.pB;
   pD = params.price.pD;
   F_0= params.price.F_0;
   
   KcB = params.gain.KcB;
   KcD = params.gain.KcD;
   MDs = params.gain.MDs;
   MBs = params.gain.MBs;
   Ds  = params.gain.Ds;
   Bs  = params.gain.Bs;
   
   C = params.colloc.C;
   D = params.colloc.D;
   h = params.colloc.h;
   
   delta_time = params.weight.delta_time;
   Qmax = params.Qmax;
   
   global nx nu nk d ns;

   for k=0:nk-1
        % New NLP variable for the control
        %Uk  = MX.sym(['U_' num2str(k)], nu);
        Uk     = MX.sym(['U_' num2str((iter-1)*nk+k)], nu);
        w      = {w{:}, Uk};
        lbw    = [lbw; u_min];
        ubw    = [ubw; u_max];
        indexU = (iter-1)*nk + (k+1);
        w0     = [w0;  u(:,indexU)];
        
      
        % State at collocation points
        Xkj   = {};
        SumX1 = 0;
        for j=1:d
            Xkj{j} = MX.sym(['X_' num2str((iter-1)*nk+k) '_' num2str(j)], nx);
            w      = {w{:}, Xkj{j}};
            lbw    = [lbw; x_min];
            ubw    = [ubw; x_max];
            w0     = [w0; x(iter+1,:)'];
            count  = count + 1;
        end

        % Loop over collocation points
        Xk_end = D(1)*Xk; 

        for j=1:d
           % Expression for the state derivative at the collocation point
           xp = C(1,j+1)*Xk;
           for r=1:d
               xp = xp + C(r+1,j+1)*Xkj{r};
           end

           % Append collocation equations
           fj  = f(Xkj{j},Uk);
           g   = {g{:}, h*fj - xp};
           lbg = [lbg; zeros(nx,1)];
           ubg = [ubg; zeros(nx,1)];

           % Add contribution to the end state
           Xk_end = Xk_end + D(j+1)*Xkj{j};
           
        end    

        % New NLP variable for state at end of interval
        Xk  = MX.sym(['X_' num2str((iter-1)*nk+k)], nx);
        w   = {w{:}, Xk};
        lbw = [lbw; x_min];
        ubw = [ubw; x_max];
        w0  = [w0; x(iter+1,:)'];
        count  = count + 1;

        % Add equality constraint
        g   = {g{:}, (Xk_end-Xk)};
        lbg = [lbg; zeros(nx,1)];
        ubg = [ubg; zeros(nx,1)];
               
        Jecon  = (pf*F_0 + pV*Uk(2) - pB*Uk(5) - pD*Uk(4)) * delta_time;
        
        alpha  = 1;
        beta   = 1;
        gamma  = 1;

        J = beta*Jecon;
        
    end
end


