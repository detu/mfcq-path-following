function [solver,par] = buildNLP5(f,par)
import casadi.*

GOR_case = par.GOR_case;
nx = par.nx;
nz = par.nz;
nu = par.nu;
dx0 = par.dx0;
z0 = par.z0;
[dx0,z0,u0,lbx,lbz,lbu,ubx,ubz,ubu] = GasLift_Initialization_bounds(par);
ubu = [10;10];
%% Direct Collocation

% Degree of interpolating polynomial
d = 3;

% Get collocation points
tau_root = [0, collocation_points(d, 'radau')];

% Coefficients of the collocation equation
C = zeros(d+1,d+1);

% Coefficients of the continuity equation
D = zeros(d+1, 1);

% Coefficients of the quadrature function
B = zeros(d+1, 1);

% Construct polynomial basis
for j=1:d+1
    % Construct Lagrange polynomials to get the polynomial basis at the collocation point
    coeff = 1;
    for r=1:d+1
        if r ~= j
            coeff = conv(coeff, [1, -tau_root(r)]);
            coeff = coeff / (tau_root(j)-tau_root(r));
        end
    end
    % Evaluate the polynomial at the final time to get the coefficients of the continuity equation
    D(j) = polyval(coeff, 1.0);
    
    % Evaluate the time derivative of the polynomial at all collocation points to get the coefficients of the continuity equation
    pder = polyder(coeff);
    for r=1:d+1
        C(j,r) = polyval(pder, tau_root(r));
    end
    
    % Evaluate the integral of the polynomial to get the coefficients of the quadrature function
    pint = polyint(coeff);
    B(j) = polyval(pint, 1.0);
end


% Preparing for scenario tree
M = size(GOR_case);
nR = 1;
nS = M(2)^nR;

a = 1;
ac = 1;
for i = 1:par.N
    for j = 1:nS
        
        nonant(i,j) = a;
        
        if i == 1
            a=1;
            
        else if i<=nR
                ac = ac+1;
                if ac>M
                    a = a+1;
                    ac = 1;
                end
            else
                a =NaN;
            end
        end
        
        
    end
end
nonant = nonant';

%% build NLP solver
% empty nlp
w = {};
w0 = [];
lbw = [];
ubw = [];
J = 0;

g = {};
lbg = [];
ubg = [];

ws = MX.sym('ws',nS);

% add steady-state variable
% incorporate steady-state optimization
[w0s,lbws,ubws,lbgs,ubgs,gs,wss] = gasLiftSSVariables(par);
g      = {g{:}, gs};
lbg    = [lbg; lbgs];
ubg    = [ubg; ubgs];
w      = {w{:}, wss};
w0     = [w0; w0s];
lbw    = [lbw; lbws];
ubw    = [ubw; ubws];

for js = 1:nS
    %     GOR_val = GOR_case(:,js);
    % initial conditions for each scenario
    X0 = MX.sym(['X0' '_' num2str(js)],nx);
    w = {w{:}, X0};
    lbw = [lbw; lbx];
    ubw = [ubw; ubx];
    w0 = [w0; dx0];
    
    % Formulate NLP
    Xk = X0;
    Xkj = {};
    Zkj = {};
    
    
    
    X0_par = MX.sym(['X0_par_' num2str(js)],nx);
    
    
    
    switch js
        case 1
            
            U0_1     = MX.sym(['U0_' num2str(js)],nu);
            Uk_prev = U0_1;
            
            X0_par_1 = MX.sym(['X0_par_' num2str(js)],nx);
            g   = {g{:},X0-X0_par_1};
            lbg = [lbg;zeros(nx,1)];
            ubg = [ubg;zeros(nx,1)];
            
            GOR_val1 = MX.sym('GOR_val1',2);
            
        case 2
            
            U0_2     = MX.sym(['U0_' num2str(js)],nu);
            Uk_prev = U0_2;
            
            X0_par_2 = MX.sym(['X0_par_' num2str(js)],nx);
            g   = {g{:},X0-X0_par_2};
            lbg = [lbg;zeros(nx,1)];
            ubg = [ubg;zeros(nx,1)];
            
            GOR_val2 = MX.sym('GOR_val2',2);
        case 3
            
            U0_3     = MX.sym(['U0_' num2str(js)],nu);
            Uk_prev = U0_3;
            
            X0_par_3 = MX.sym(['X0_par_' num2str(js)],nx);
            g   = {g{:},X0-X0_par_3};
            lbg = [lbg;zeros(nx,1)];
            ubg = [ubg;zeros(nx,1)];
            
            GOR_val3 = MX.sym('GOR_val3',2);
        case 4
            
            U0_4     = MX.sym(['U0_' num2str(js)],nu);
            Uk_prev = U0_4;
            
            X0_par_4 = MX.sym(['X0_par_' num2str(js)],nx);
            g   = {g{:},X0-X0_par_4};
            lbg = [lbg;zeros(nx,1)];
            ubg = [ubg;zeros(nx,1)];
            
            GOR_val4 = MX.sym('GOR_val4',2);
            
        case 5
            
            U0_5     = MX.sym(['U0_' num2str(js)],nu);
            Uk_prev = U0_5;
            
            X0_par_5 = MX.sym(['X0_par_' num2str(js)],nx);
            g   = {g{:},X0-X0_par_5};
            lbg = [lbg;zeros(nx,1)];
            ubg = [ubg;zeros(nx,1)];
            
            GOR_val5 = MX.sym('GOR_val5',2);
    end
    
    for k = 0:par.N-1
        
        Uk  = MX.sym(['U_' num2str(k) '_' num2str(js)],nu);
        
        w   = {w{:},Uk};
        lbw = [lbw;lbu];
        ubw = [ubw;ubu];
        w0  = [w0;u0];
        
        
        Xkj = {};
        Zkj = {};
        
        for j = 1:d
            Xkj{j} = MX.sym(['X_' num2str(k) '_' num2str(j) '_' num2str(js)],nx);
            Zkj{j} = MX.sym(['Z_' num2str(k) '_' num2str(j) '_' num2str(js)],nz);
            
            w   = {w{:},Xkj{j},Zkj{j}};
            lbw = [lbw;lbx;lbz];
            ubw = [ubw;ubx;ubz];
            w0  = [w0; dx0;z0 ];
        end
        
        % Loop over collocation points
        Xk_end  = D(1)*Xk;
        
        for j = 1:d
            % Expression for the state derivative at the collocation point
            xp  = C(1,j+1)*Xk;  % helper state
            for r = 1:d
                xp = xp + C(r+1,j+1)*Xkj{r};
            end
            
            switch js
                case 1
                    [fj,zj,qj] =  f(Xkj{j},Zkj{j},vertcat(Uk,GOR_val1));
                case 2
                    [fj,zj,qj] =  f(Xkj{j},Zkj{j},vertcat(Uk,GOR_val2));
                case 3
                    [fj,zj,qj] =  f(Xkj{j},Zkj{j},vertcat(Uk,GOR_val3));
                case 4
                    [fj,zj,qj] =  f(Xkj{j},Zkj{j},vertcat(Uk,GOR_val4));
                case 5
                    [fj,zj,qj] =  f(Xkj{j},Zkj{j},vertcat(Uk,GOR_val5));
            end
            
            g   = {g{:},par.tf*fj-xp,zj};  % dynamics and algebraic constraints
            lbg = [lbg;zeros(nx,1);zeros(nz,1)];
            ubg = [ubg;zeros(nx,1);zeros(nz,1)];
            
            % Gas capacity constraints on all the collocation points
            g = {g{:},sum(Zkj{j}(17:18))};
            lbg = [lbg;0];
            ubg = [ubg;par.QgMax];
            
            % Add contribution to the end states
            Xk_end = Xk_end + D(j+1)*Xkj{j};
            
            %J = J + ws(js)*(B(j+1)*qj*par.tf ) + 10*(sum(Uk-Uk_prev).^2);
            J = J + ws(js)*(B(j+1)*qj*par.tf ) + 10*(sum(Uk-wss{3}).^2);
            
            
        end
%         g = {g{:},Uk_prev-Uk};
%         lbg = [lbg;-1;-1];
%         ubg = [ubg;1;1];
%         
%         Uk_prev = MX.sym(['Uprev_' num2str(k+1) '_' num2str(js)],nu);
%         Uk_prev = Uk;
        
        % New NLP variable for state at end of interval
        Xk = MX.sym(['X_' num2str(k+1) '_' num2str(js)], nx);
        
        w = {w{:},Xk};
        lbw = [lbw;lbx];
        ubw = [ubw;ubx];
        w0 = [w0; dx0];
        
        % Shooting Gap constraint
        g = {g{:},Xk_end-Xk};
        lbg = [lbg;zeros(nx,1)];
        ubg = [ubg;zeros(nx,1)];
        
        g = {g{:},sum(Uk)};
        lbg = [lbg;0];
        ubg = [ubg;100.*par.qGLMax];
        if k < nR
            u_ant{js,k+1} = MX.sym(['u_ant_' num2str(k) '_' num2str(js)],nu);
            u_ant{js,k+1} = Uk;
        end
        
    end
    
end

% Add Non-anticipativity constraints
for k = 1:nR
    for js = 1:nS-1
        if ~isnan(nonant(js,k))
            if nonant(js) == nonant(js+1)
                g = {g{:},(u_ant{js,k} - u_ant{js+1,k})};
                lbg = [lbg;zeros(nu,1)];
                ubg = [ubg;zeros(nu,1)];
            end
        end
    end
end

% create and solve NLP solver

opts = struct('warn_initial_bounds',false, ...
    'print_time',false, ...
    'ipopt',struct('print_level',5) ...
    );


nlp = struct('x',vertcat(w{:}),'p',vertcat(GOR_val1,GOR_val2,GOR_val3,GOR_val4,GOR_val5,...
    U0_1,U0_2,U0_3,U0_4,U0_5,X0_par_1,X0_par_2,X0_par_3,X0_par_4,X0_par_5,ws),'f',J,'g',vertcat(g{:}));

solver = nlpsol('solver','ipopt',nlp,opts);

par.w0 = w0;
par.lbw = lbw;
par.ubw = ubw;
par.lbg = lbg;
par.ubg = ubg;
par.nS = nS;