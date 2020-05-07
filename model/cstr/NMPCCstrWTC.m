function NMPCCstr
%NMPCCSTR Summary of this function goes here
% 
% NMPC controller for CSTR 
%
% [OUTPUTARGS] = NMPCCSTR(INPUTARGS) Explain usage here
% 
% Examples: 
% 
% Provide sample usage code here
% 
% See also: List related files here

% $Author: suwartad $	$Date: 2017/10/07 15:37:22 $	$Revision: 0.1 $
% Copyright: Process Control Group - NTNU Trondheim 2017

% TO DO:
% - add regularization term for state
% - compute weights based on Gershgorin bound criteria

import casadi.* 
format long;

global N k pQ pC;
% number of mpc iteration
mpciterations = 100;
% number of prediction horizon
N             = 50;
%N             = 1;
% sampling time
T             = 1;  % [minute]
% initial control and state
u0            = 20;
u0            = repmat(u0,1,N);
% get initial measurement (states) at time T = 0.
tmeasure      = 0.0;
xmeasure      = [0.9;0.6];


% new code dedicated for CSTR
[~, xmeasureAll, uAll, obj, optRes, params, runtime] = iNmpcCstr(@optProblem, @system, mpciterations, N, T, tmeasure, xmeasure, u0);


keyboard;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Definition of the NMPC functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = system(t, x, u, T)
    
    global uc;
    uc = u;
    %[~,x_out] = ode15s('cola_lv_cstr',[t t+T], x);
    [~,x_out] = ode15s('cstr_odeWTC',[t t+T], x);
    lengthx   = size(x_out); 
    y         = x_out(lengthx(1),:)'; 
    
end

function [J,g,w0,w,lbg,ubg,lbw,ubw,params] = optProblem(x, u, N, x0_measure)   %add prediction horizon 
    import casadi.*
    
    % invoke the model
    [~,state,xdot,inputs] = cstrFuncWTC();
    f = Function('f',{state,inputs}, {xdot});
    
    % bound constraints
    Q_max = 20;
    Q_min = 10;
    
    
    % State bounds and initial guess
    x_min =  zeros(2,1);  
    x_max =  ones(2,1);
    
    % state constraint on Cb
    %x_min(2)  = 0.45;
    %x_min(2)  = 0.49;
    %x_max(2)  = 1;
    
    % Control bounds
    u_min = Q_min;
    u_max = Q_max;
    
    % compact bound variable for ease of function invocation 
    params.bound.x_min = x_min;
    params.bound.x_max = x_max;
    params.bound.u_min = u_min;
    params.bound.u_max = u_max;
    
    % Construct objective function
    load CstrXinit.mat
    xf    = Xinit(1:2);
    u_opt = Xinit(3);
    
    % dimensions
    global nx nu nk d tf ns;
    nx = 2;   
    nu = 1;    
    nk = 1;
    %tf = 1;   % in [minutes]
    tf = .1;
    %tf = 3;    % 3 minutes
    h  = tf/nk;
    ns = 0;
    
    % compact model variable
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
    
    % incorporate steady-state optimization
    u1      = MX.sym('u1');     % u  
    x1      = MX.sym('x1');     % x1 = Ca
    x2      = MX.sym('x2');     % x2 = Cb
    V       = 10;
    Caf     = 1;
    %epsilon = MX.sym('epsilon');
    xs     = [x1;x2;u1]; 
    global k;
    x1dot  = (u1/V)*(Caf - x1) - k*x1;
    x2dot  = (u1/V)*-x2 + k*x1;
    ceq    = [x1dot;x2dot];
    g      = {g{:}, ceq};
    lbg    = [lbg; 0; 0];
    ubg    = [ubg; 0; 0];
    w      = {w{:}, xs};
    w0     = [w0; x0_measure; u(1)];
    lbx    = [0;0;10];
    ubx    = [1;1;20];
    lbw    = [lbw; lbx];
    ubw    = [ubw; ubx];
    
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

%     %load Qmax.mat;
%     load QmaxCstr.mat;
%     params.Qmax = Qmax;
    
    % HERE SHOULD BE LOOP N-TIMES ACCORDING TO THE NUMBER OF PREDICTION HORIZON
    count  = 2; % counter for state variable as initial guess
    ssoftc = 0;
    for i=1:N
        [J,g,w0,w,lbg,ubg,lbw,ubw,Xk,params,count,ssoftc] = iterateOnPredictionHorizon(Xk, w, w0, lbw, ubw, lbg, ubg, g, J, params, i, count, ssoftc, xs);
    end

    
end


function [J,g,w0,w,lbg,ubg,lbw,ubw,Xk,params,count,ssoftc] = iterateOnPredictionHorizon(Xk, w, w0, lbw, ubw, lbg, ubg, g, J, params, iter, count, ssoftc, xs)

   import casadi.*
   global N pC pQ;
   % extract compact variables
   x_min = params.bound.x_min;
   x_max = params.bound.x_max;
   u_min = params.bound.u_min;
   u_max = params.bound.u_max;
   
   f  = params.model.f;
   xdot_val_rf_ss = params.model.xdot_val_rf_ss;
   x = params.model.x;
   u = params.model.u;
   u_opt = params.model.u_opt;
  
   
   C = params.colloc.C;
   D = params.colloc.D;
   h = params.colloc.h;
   
   delta_time = params.weight.delta_time;
   %Qmax = params.Qmax;
   
   global nx nu nk d ns xf uf check;

   for k=0:nk-1
        % New NLP variable for the control
        %Uk  = MX.sym(['U_' num2str(k)], nu);
        Uk     = MX.sym(['U_' num2str((iter-1)*nk+k)], nu);
        w      = {w{:}, Uk};
        lbw    = [lbw; u_min];
        ubw    = [ubw; u_max];
        indexU = (iter-1)*nk + (k+1);
        w0     = [w0;  u(:,indexU)];
        
%         %Jcontrol   = (Qmax(nx+1:nx+nu,1).*(Uk - u_opt))' * (Uk - u_opt);
%         if check == 0
%             Jcontrol   = 2*(Uk - xs(3)) * (Uk - xs(3));
%         else
%             Jcontrol   = 2*(Uk - uf) * (Uk - uf);
%         end
        
        Jcontrol   = 2*(Uk - xs(3)) * (Uk - xs(3));
      
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
               
        %Jecon  = -Uk(1)*(2*Xk(2) - 0.5)* delta_time;
        %Jecon  = -Uk(1)* delta_time;
        
        Jecon  = -Uk(1)*(pC*Xk(2) - pQ)* delta_time;
%         Jstate =(Qmax(1:nx,1).*(Xk - xdot_val_rf_ss))' * (Xk - xdot_val_rf_ss) * delta_time;
%         
%         alpha  = 1;
%         beta   = 1;
%         gamma  = 1;

        
%         % compute rotated cost function
%         fm  = f(Xk,Uk);
%         % load Lagrange multipliers from steady-state optimization
%         load LamdaCstrDist.mat; % lamda
%         Jmodel = lamda'*fm;
        
        %J = J + alpha*Jcontrol + gamma*Jstate + beta*Jecon;
        %J = J + Jecon;
        J = J + Jecon + Jcontrol;
    end
end