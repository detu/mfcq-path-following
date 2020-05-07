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
load Xinit3075.mat
u0            = Xinit3075(85:89);
u0            = repmat(u0,1,N);
% get initial measurement (states) at time T = 0.
tmeasure      = 0.0;
xmeasure      = Xinit3075(1:84);

load Xopt3175Horizon90.mat
xGuess = Xopt;

% new iNmpcDual that collect dual variable as well 
[~, xmeasureAll, uAll, obj, primalRes, dualRes, params, runtime] = iNmpcDual(@optProblem, @system, mpciterations, N, T, tmeasure, xmeasure, u0, xGuess);

save results.mat xmeasureAll uAll dualRes runtime;

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

function [J,g,w0,w,lbg,ubg,lbw,ubw,params] = optProblem(x, u, N, x0_measure)    %add prediction horizon 
    import casadi.*
    
    %xc = 10;
    xc = 1;
    %scaling for bottom concentration
    %x(:,1) = xc/10;
    
    % the model
    NT = 41;
    global Uf;
    Uf = 0.31;           % Feeding rate F_0
    
    % invoke the model
    [~,state,xdot,inputs] = DistColACstr(Uf);
    %[~,state,xdot,inputs] = DistColACstrScaled(Uf);  % with scaled equations
    %f = Function('f',{state,inputs}, {xdot});
    
    % objective function value
    obj = computeObjectiveFunction(state,inputs);
 
    f = Function('f',{state,inputs}, {xdot,obj});
    
    % bound constraints
    VB_max = 4.008;
    xB_max = xc/10;    %scaled bottom concentration
    
    
    % State bounds and initial guess
    x_min =  zeros(84,1);  % try without epsilon here, later put epsilon
    x_max =  ones(84,1);
    
    %x_min(84) = 0.3;
    x_max(1)  = xB_max; % scaled bottom concentration
    x_max(84) = 0.75;

    
    % Control bounds
    u_min = [0.1; 0.1; 0.1; 0.1; 0.1];
    u_max = [10; VB_max; 10; 1.0; 1.0];
    
    % soft constraint bounds
    sc_min = 0;
    sc_max = 1e6;
    
    % compact bound variable for ease of function invocation 
    params.bound.x_min  = x_min;
    params.bound.x_max  = x_max;
    params.bound.u_min  = u_min;
    params.bound.u_max  = u_max;
    params.bound.sc_min = sc_min;
    params.bound.sc_max = sc_max;
    
    % Construct objective function
    load CstrDistXinit.mat;
    xf    = Xinit(1:84);
    xf(1) = xc*xf(1);          % scaled bottom concentration
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
    [B,C,D,d] = collocationSetup();
    
    % compact collocation variable
    params.colloc.C = C;
    params.colloc.D = D;
    params.colloc.h = h;
    params.colloc.B = B;

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
    
    %w0  = [w0; x(1,1:nx)'];
    w0  = [w0; x0_measure];
    
    g   = {g{:}, X0 - x0_measure};  % USE MEASUREMENT HERE !
    lbg = [lbg; zeros(nx,1)];
    ubg = [ubg; zeros(nx,1)];

    % formulate the NLP
    Xk = X0;

    load Qmax.mat;
    params.Qmax = Qmax;
   
    
    % HERE SHOULD BE LOOP N-TIMES ACCORDING TO THE NUMBER OF PREDICTION HORIZON
    %count  = 2; % counter for state variable as initial guess
    count = 1;
    ssoftc = 0;
    for i=1:N
        [J,g,w0,w,lbg,ubg,lbw,ubw,Xk,params,count,ssoftc] = iterateOnPredictionHorizon(Xk, w, w0, lbw, ubw, lbg, ubg, g, J, params, i, count,ssoftc);
    end

    
end


function [J,g,w0,w,lbg,ubg,lbw,ubw,Xk,params,count,ssoftc] = iterateOnPredictionHorizon(Xk, w, w0, lbw, ubw, lbg, ubg, g, J, params, iter, count, ssoftc)

   import casadi.*
   global N;
   % extract compact variables
   x_min  = params.bound.x_min;
   x_max  = params.bound.x_max;
   u_min  = params.bound.u_min;
   u_max  = params.bound.u_max;
   sc_min = params.bound.sc_min;
   sc_max = params.bound.sc_max;
   
   
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
   
   C = params.colloc.C;
   D = params.colloc.D;
   h = params.colloc.h;
   B = params.colloc.B;
   
   delta_time = params.weight.delta_time;
   Qmax = params.Qmax;
   
   global nx nu nk d tf ns;
   

% with collocation method
   for k=0:nk-1
        % New NLP variable for the control
        Uk     = MX.sym(['U_' num2str((iter-1)*nk+k)], nu);
        w      = {w{:}, Uk};
        lbw    = [lbw; u_min];
        ubw    = [ubw; u_max];
        indexU = (iter-1)*nk + (k+1);
        w0     = [w0;  u(:,indexU)];
        
        Jcontrol   = (Qmax(nx+1:nx+nu,1).*(Uk - u_opt))' * (Uk - u_opt);
      
        % State at collocation points
        Xkj   = {};
        Skj   = {};
        Ekj   = {};
        Jcoll = 0;
        sumE  = {};
        for j=1:d
            Xkj{j} = MX.sym(['X_' num2str((iter-1)*nk+k) '_' num2str(j)], nx);
            w      = {w{:}, Xkj{j}};
            lbw    = [lbw; x_min];
            ubw    = [ubw; x_max];

            w0     = [w0; x(:,count)];
            count  = count + 1;
                        
            Jcoll = Jcoll + (Qmax(1:nx,1).*(Xkj{j} - xdot_val_rf_ss))' * (Xkj{j} - xdot_val_rf_ss) * delta_time;

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
           %fj  = f(Xkj{j},Uk);
           [fj, qj] = f(Xkj{j},Uk);
           g   = {g{:}, h*fj - xp};
           lbg = [lbg; zeros(nx,1)];
           ubg = [ubg; zeros(nx,1)];

           % Add contribution to the end state
           Xk_end = Xk_end + D(j+1)*Xkj{j};
           
           %J = J + B(j+1)*qj*h;
        end    

        % New NLP variable for state at end of interval
        Xk  = MX.sym(['X_' num2str((iter-1)*nk+k)], nx);
        w   = {w{:}, Xk};
        
        lbw = [lbw; x_min];
        ubw = [ubw; x_max];
        w0   = [w0; x(:,count)];
        count  = count + 1;

        % Add equality constraint
        g   = {g{:}, (Xk_end-Xk)};
        lbg = [lbg; zeros(nx,1)];
        ubg = [ubg; zeros(nx,1)];
       
               
        Jecon  = (pf*F_0 + pV*Uk(2) - pB*Uk(5) - pD*Uk(4)) * delta_time;
        Jstate =(Qmax(1:nx,1).*(Xk - xdot_val_rf_ss))' * (Xk - xdot_val_rf_ss) * delta_time;
        
        alpha  = 0;
        beta   = 1;
        gamma  = 1;

        J = J + alpha*Jcontrol + gamma*Jstate + beta*Jecon + Jcoll;
		%J = J + alpha*Jcontrol + gamma*Jstate + beta*Jecon;
        
    end

end


