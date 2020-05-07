function [prob] = cstrMfcq(p)
%CSTRMFCQ Summary of this function goes here
% 
% [OUTPUTARGS] = CSTRMFCQ(INPUTARGS) Explain usage here
% 
% Examples: 
% 
% Provide sample usage code here
% 
% See also: List related files here

% $Author: suwartad $	$Date: 2017/10/07 21:50:54 $	$Revision: 0.1 $
% Copyright: Process Control Group - NTNU Trondheim 2017

import casadi.*

prob = struct('neq',{0},'niq',{0},'cin',{0},'ceq',{0},'dp_in',{0},'dp_eq',{0},'hess',{0},'lxp',{0},'x',0,'name',0);
prob.neq  = 2000;         % HARD-CODE !    % number of equality constraint
prob.niq  = 0;            % number of inequality constraint
prob.name = 'CSTR';
prob.x    = zeros(2,1);

prob.obj  = (@(x,y,p,N)(objective(x,y,p,N)));

end

function [f,g,H,Lxp,cst,J,cp,Jeq,dpe,Hobj,Hc] = objective(x,y,p,N)
    import casadi.*
   
    nPrimal = numel(x);
    nDual   = numel(y.lam_g);
    nParam  = numel(p);
    % model parameters
    [~,state,xdot,inputs] = cstrFunc();
    sf = Function('sf',{state,inputs}, {xdot});
    
    % dimensions
    global nx nu nk d tf ns;
    h  = tf/nk;
    
    % preparing collocation matrices
    [~,C,D,d] = collocationSetup();
    

    % NLP variable vector 
    V    = {};      % decision variables contain both control and state variables
    obj  = 0;       % objective function
    cons = {};      % nonlinear constraint
    
    delta_time = 1;
    alpha = 1;
    beta  = 1;
    gamma = 1;

    % Initial states and controls 
    load CstrXinit.mat
    xf    = Xinit(1:2);
    u_opt = Xinit(3);
    
    % compact model variable
    params.model.sf  = sf;
    params.model.xdot_val_rf_ss = xf;
    params.model.x  = x;
    params.model.u_opt = u_opt;
    
    % compact collocation variable
    params.colloc.C = C;
    params.colloc.D = D;
    params.colloc.h = h;
    
    % compact weight variable
    params.weight.delta_time = delta_time;
    params.weight.alpha      = alpha;
    params.weight.beta       = beta;
    params.weight.gamma      = gamma;
    
    % "Lift" initial conditions
    X0      = MX.sym('X0', nx);
    V       = {V{:}, X0};
    cons    = {cons{:}, X0 - x(1:nx,1)};
    cons_x0 = X0 - x(1:nx,1);
    
    % formulate the NLP
    Xk = X0;

    load QmaxCstr.mat;
    params.Qmax = Qmax;
    ssoftc = 0;
    for i=1:N
        [obj,cons,V,Xk,params,ssoftc] = iterateOnPredictionHorizon(Xk, V, cons, obj, params, i,ssoftc);
    end
    
    V = vertcat(V{:});
    % Concatenate constraints
    cons  = vertcat(cons{:});
    
    % objective function and constraint functions
    f = Function('f', {V}, {obj}, char('V'), char('objective'));
    c = Function('c', {V}, {cons}, char('V'), char('constraint'));
    cx0 = Function('cx0', {X0}, {cons_x0}, char('X0'), char('constraint'));
    
    % construct Lagrangian
    lag_expr = obj + y.lam_g'*cons;  
    
    g    = f.gradient();
    lagr = Function('lagr', {V}, {lag_expr}, char('V'), char('lag_expr'));
    H    = Function(lagr.hessian('V','lag_expr'));
    Hobj = f.hessian('V','objective');
    J    = c.jacobian('V','constraint');
    Jp   = cx0.jacobian('X0','constraint');
   
    global flagDt;
    if (flagDt == 1)
        %tic;
        for i=1:nDual
            % return in cell array
            ci        = Function('ci', {V}, {cons(i,1)}, char('V'), char('ci_expr'));
            Hci       = Function(ci.hessian('V','ci_expr'));
            Hc{:,:,i} = sparse(Hci(x));
        end
        %toc;
        save HcCstr.mat Hc;
    else
        Hc = [];
    end
    
    f                        = f(x);
    g                        = g(x);
    H                        = H(x);
    Lxp                      = H(1:nPrimal,1:nParam);
    J                        = J(x);
    Jtemp                    = zeros(nDual,nParam);
    cp                       = Jp(x(1:nParam));
    Jtemp(1:nParam,1:nParam) = full(cp);
    cp                       = sparse(Jtemp);
    cst                      = c(x);

    
    % Evaluation of objective function used for Greshgorin bound
    Hobj = Hobj(x);
    Hobj = sparse(Hobj);

    f   = full(f);
    g   = sparse(g);    
    H   = sparse(H);
    Lxp = sparse(Lxp);
    J   = sparse(J);
    cp  = sparse(cp);
    cst = full(cst);

    % Equality constraint
    Jeq = J;
    dpe = cp;
    
    
end

function [obj,cons,V,Xk,params,ssoftc] = iterateOnPredictionHorizon(Xk, V, cons, obj, params, iter, ssoftc)

   import casadi.*
   
   % extract compact variables

   sf  = params.model.sf;
   xdot_val_rf_ss = params.model.xdot_val_rf_ss;
   u_opt = params.model.u_opt;
   

   C = params.colloc.C;
   D = params.colloc.D;
   h = params.colloc.h;
   
   delta_time = params.weight.delta_time;
   Qmax = params.Qmax;
   
   global nx nu nk d;

   for k=0:nk-1
        % New NLP variable for the control
        %Uk  = MX.sym(['U_' num2str(k)], nu);
        Uk  = MX.sym(['U_' num2str((iter-1)*nk+k)], nu);
        V   = {V{:}, Uk};
        
        Jcontrol   = (Qmax(nx+1:nx+nu,1).*(Uk - u_opt))' * (Uk - u_opt);

        % State at collocation points
        Xkj = {};
        for j=1:d
            %Xkj{j} = MX.sym(['X_' num2str(k) '_' num2str(j)], nx);
            Xkj{j} = MX.sym(['X_' num2str((iter-1)*nk+k) '_' num2str(j)], nx);
            V      = {V{:}, Xkj{j}};
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
           fj   = sf(Xkj{j},Uk);
           cons = {cons{:}, h*fj - xp};

           % Add contribution to the end state
           Xk_end = Xk_end + D(j+1)*Xkj{j};

        end    

        % New NLP variable for state at end of interval
        Xk  = MX.sym(['X_' num2str((iter-1)*nk+k)], nx);
        V   = {V{:}, Xk};

        % Add equality constraint
        cons= {cons{:}, Xk_end-Xk};

        %Jecon  = -Uk(1)*(2*Xk(2) - 0.5)* delta_time;
        Jecon  = -Uk(1)* delta_time;
        Jstate =(Qmax(1:nx,1).*(Xk - xdot_val_rf_ss))' * (Xk - xdot_val_rf_ss) * delta_time;
        
        alpha  = 1;
        beta   = 1;
        gamma  = 1;
        obj    = obj + alpha*Jcontrol + gamma*Jstate + beta*Jecon;

    end
end