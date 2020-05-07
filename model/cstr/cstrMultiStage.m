function [prob] = cstrMultiStage(p)
%CSTRMULTISTAGE Summary of this function goes here
% 
% [OUTPUTARGS] = CSTRMULTISTAGE(INPUTARGS) Explain usage here
% 
% Examples: 
% 
% Provide sample usage code here
% 
% See also: List related files here

% $Author: suwartad $	$Date: 2018/05/04 08:37:29 $	$Revision: 0.1 $
% Copyright: Process Control Group - NTNU Trondheim 2018

import casadi.*

prob = struct('neq',{0},'niq',{0},'cin',{0},'ceq',{0},'dp_in',{0},'dp_eq',{0},'hess',{0},'lxp',{0},'x',0,'name',0);
prob.neq  = 2000;         % HARD-CODE !    % number of equality constraint
prob.niq  = 0;            % number of inequality constraint
prob.name = 'CSTR';
prob.x    = zeros(2,1);

prob.obj  = (@(x,y,p,N)(objective(x,y,p,N)));

end

function [f,g,H,Lxp,cst,J,cp,Jeq,dpe,Hobj] = objective(x,y,p,N)
    import casadi.*
   
    nPrimal = numel(x);
    nDual   = numel(y.lam_g);
    nParam  = numel(p);
    % model parameters
    [~,state,xdot,inputs] = cstrFunc();
    sf = Function('sf',{state,inputs}, {xdot});
    
    % dimensions
    global nx nk d tf k scenario;
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
    
    % compact model variable
    params.model.sf = sf;
    params.model.x  = x;
    
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
    V0      = X0;
    cons    = {cons{:}, X0 - x(1:nx,1)};
    cons_x0 = X0 - x(1:nx,1);
    
    k = 1.2;
    if scenario == 1
        % incorporate steady-state optimization
        u1_1   = MX.sym('u1_1');     % u
        x1_1   = MX.sym('x1_1');     % x1 = Ca
        x2_1   = MX.sym('x2_1');     % x2 = Cb
        xs     = [x1_1;x2_1;u1_1];
        x1dot  = (u1_1/10)*(1 - x1_1) - k*x1_1;
        x2dot  = (u1_1/10)*-x2_1 + k*x1_1;
    end
    
    if scenario == 2
        % incorporate steady-state optimization
        u1_2   = MX.sym('u1_2');     % u
        x1_2   = MX.sym('x1_2');     % x1 = Ca
        x2_2   = MX.sym('x2_2');     % x2 = Cb
        xs     = [x1_2;x2_2;u1_2];
        x1dot  = (u1_2/10)*(1 - x1_2) - k*x1_2;
        x2dot  = (u1_2/10)*-x2_2 + k*x1_2;
    end
    
    if scenario == 3
        % incorporate steady-state optimization
        u1_3   = MX.sym('u1_3');     % u
        x1_3   = MX.sym('x1_3');     % x1 = Ca
        x2_3   = MX.sym('x2_3');     % x2 = Cb
        xs     = [x1_3;x2_3;u1_3];
        x1dot  = (u1_3/10)*(1 - x1_3) - k*x1_3;
        x2dot  = (u1_3/10)*-x2_3 + k*x1_3;
    end
    
    if scenario == 4
        % incorporate steady-state optimization
        u1_4   = MX.sym('u1_4');     % u
        x1_4   = MX.sym('x1_4');     % x1 = Ca
        x2_4   = MX.sym('x2_4');     % x2 = Cb
        xs     = [x1_4;x2_4;u1_4];
        x1dot  = (u1_4/10)*(1 - x1_4) - k*x1_4;
        x2dot  = (u1_4/10)*-x2_4 + k*x1_4;
    end
    

%     % incorporate steady-state optimization
%     u1     = MX.sym('u1');     % u
%     x1     = MX.sym('x1');     % x1 = Ca
%     x2     = MX.sym('x2');     % x2 = Cb
%     xs     = [x1;x2;u1];
%     x1dot  = (u1/10)*(1 - x1) - k*x1;
%     x2dot  = (u1/10)*-x2 + k*x1;
%     
    % constraint of the artificial steady-state 
    ceq    = [x1dot;x2dot];
    cons   = {cons{:}, ceq};
    V      = {V{:}, xs};
    
%     % "Lift" initial conditions
%     X0      = MX.sym('X0', nx);
%     V0      = X0;
%     V       = {V{:}, X0};
%     cons    = {cons{:}, X0 - x(1:nx,1)};
%     cons_x0 = X0 - x(1:nx,1);
    
    % formulate the NLP
    Xk = X0;

    ssoftc = 0;
    for i=1:N
        [obj,cons,V,Xk,params,ssoftc] = iterateOnPredictionHorizon(Xk, V, cons, obj, params, i, ssoftc, xs);
    end
    
    V     = vertcat(V{:});
    % Concatenate constraints
    cons  = vertcat(cons{:});
    
    % objective function and constraint functions
    f = Function('f', {V}, {obj}, char('V'), char('objective'));
    c = Function('c', {V}, {cons}, char('V'), char('constraint'));
    cx0 = Function('cx0', {X0}, {cons_x0}, char('X0'), char('constraint'));
    
    % construct Lagrangian
    lag_expr = obj + y.lam_g'*cons;  
    
    %g    = f.gradient();
    [Hj,g] = hessian(obj,V);
    Hobj   = Function('Hj',{V},{Hj});
    gj     = Function('gj',{V},{g});
%     lagr = Function('lagr', {V}, {lag_expr}, char('V'), char('lag_expr'));
%     H    = Function(lagr.hessian('V','lag_expr'));
    Hl     = hessian(lag_expr,V);
    H      = Function('Hl',{V},{Hl});
%     Hobj = f.hessian('V','objective');
%     J    = c.jacobian('V','constraint');
%     Jp   = cx0.jacobian('X0','constraint');
    J      = Function('J',{V},{jacobian(cons,V),cons});
    Jp     = Function('Jp',{V0},{jacobian(cons_x0,V0),cons_x0});
    
    f                        = f(x);
    %g                        = g(x);
    g                        = gj(x);
    H                        = H(x);
    Lxp                      = H(1:nPrimal,1:nParam);
    J                        = J(x);
    Jtemp                    = zeros(nDual,nParam);
    cp                       = Jp(x(1:nParam));
    Jtemp(1:nParam,1:nParam) = full(cp);
    cp                       = sparse(Jtemp);
    cst                      = c(x);
    
%     global flagDt;
%     if (flagDt == 1)
%         %tic;
%         for i=1:nDual
%             % return in cell array
%             %ci        = Function('ci', {V}, {cons(i,1)}, char('V'), char('ci_expr'));
%             %Hci       = Function(ci.hessian('V','ci_expr'));
%             [Hcons,~] = hessian(cons(i,1),V);
%             Hci       = Function('Hc',{V},{Hcons});
%             Hc{:,:,i} = sparse(Hci(x));
%         end
%         %toc;
%         save HcCstr.mat Hc;
%     else
%         Hc = [];
%     end

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

function [obj,cons,V,Xk,params,ssoftc] = iterateOnPredictionHorizon(Xk, V, cons, obj, params, iter, ssoftc, xs)

   import casadi.*
   
   % extract compact variables
   sf = params.model.sf;
   C  = params.colloc.C;
   D  = params.colloc.D;
   h  = params.colloc.h;
   delta_time = params.weight.delta_time;
   
   global nx nu nk d pC pQ scenario;

   for k=0:nk-1
        % New NLP variable for the control
        Uk  = MX.sym(['U_' num2str((iter-1)*nk+k)], nu);
        V   = {V{:}, Uk};
        
        % objective function with control's regularization term
        Jcontrol   = 2*(Uk - xs(3)) * (Uk - xs(3));
        %Jcontrol   = 1e-2*(Uk - xs(3)) * (Uk - xs(3));

        % State at collocation points
        Xkj = {};
        for j=1:d
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

        switch scenario 
            case 1
                Jecon = -Uk(1)*((pC+0.5*pC)*Xk(2) - (pQ+0.5*pQ))* delta_time;
            case 2
                Jecon = -Uk(1)*((pC-0.5*pC)*Xk(2) - (pQ-0.5*pQ))* delta_time;
            case 3
                Jecon = -Uk(1)*((pC+0.25*pC)*Xk(2) - (pQ+0.25*pQ))* delta_time;
            case 4
                Jecon = -Uk(1)*((pC-0.25*pC)*Xk(2) - (pQ-0.25*pQ))* delta_time;
        end

        obj = obj + Jecon + Jcontrol;
        %obj = obj + Jecon;

    end
end