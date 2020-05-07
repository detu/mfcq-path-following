
% u = 0:0.1:10;
% x1 = -(-20*exp(-0.55.*u)  -1.05.*exp(0.155.*u)+ 43)+50;
% x2 = (2.*u- 40)+50;
% plot(u,x1,u,x2)

clear
clc

% addpath ('C:\Users\dineshk\CasADi\casadi-matlabR2014b-v3.1.0-rc1')
import casadi.*


par.tf = 1;
par.N = 60;
par.nIter = 60;

%% Process model
x1 = MX.sym('x1');
x2 = MX.sym('x2');
u = MX.sym('u');
p = MX.sym('p');

tau = 20;

x1_ss = (-30*exp(-0.5.*u)  -2.*exp(0.25.*u)+ 32);
x2_ss = u+p*x1_ss;

dx1 = (x1_ss-x1)/tau;
dx2 = (x2_ss-x2)/tau;

diff = vertcat(dx1,dx2);
x_var = vertcat(x1,x2);
p_var = vertcat(u,p);

J = -x1;

ode = struct('x',x_var,'p',p_var,'ode',diff,'quad',J);
opts = struct('tf',par.tf);

% create IDAS integrator
F = integrator('F','cvodes',ode,opts);

f = Function('f',{x_var,p_var},{diff,J},{'x','p'},{'xdot','qj'});

%% Initalization

lbx = [0;0];
ubx = [100;4];
lbu = 0;
ubu = 7;
dx0 = [17.6662;3.2366];
u0 = 2;

nx = 2;
nu = 1;
%% Prepare for scenario tree generation

p_case = [0.05;0.06;0.07]; % 3 discrete realizations of p 

p_val = 0.06; % value of p used in the actual plant.

M = numel(p_case);
nR = 2;  % robust horizon 
nS = M^nR;

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


%% Optimization

% empty nlp
w = {};
w0 = [];
lbw = [];
ubw = [];
J = 0;

g = {};
lbg = [];
ubg = [];

for js = 1:nS
    % initial conditions for each scenario
    X0 = MX.sym(['X0' '_' num2str(js)],nx);
    w = {w{:}, X0};
    lbw = [lbw; lbx];
    ubw = [ubw; ubx];
    w0 = [w0; dx0];
    
    p = MX.sym('p',3);
    U0 = MX.sym('U0',1);
    X0_par = MX.sym(['X0_par''_' num2str(js)],nx);
    %
    g = {g{:},X0};
    lbg = [lbg;dx0];
    ubg = [ubg;dx0];
    
    
    % Formulate NLP
    
    Xk = X0;
    
    Xkj = {};
    Zkj = {};
    Uk_prev = U0;
    
    for k = 0:par.N-1
        
        Uk = MX.sym(['U_' num2str(k) '_' num2str(js) ],nu);
        w = {w{:},Uk};
        lbw = [lbw;lbu];
        ubw = [ubw;ubu];
        w0 = [w0;u0];
        
        % state at collocation points = s(\theta)
        Xkj = {};
        Zkj = {};
        
        for j = 1:d
            Xkj{j} = MX.sym(['X_' num2str(k) '_' num2str(j) '_' num2str(js) ],nx);
            w = {w{:},Xkj{j}};
            lbw = [lbw;lbx];
            ubw = [ubw;ubx];
            w0 = [w0; dx0];
        end
        
        % Loop over collocation points
        Xk_end = D(1)*Xk;
        
        for j = 1:d
            
            xp = C(1,j+1)*Xk;  % helper state
            
            for r = 1:d
                xp = xp + C(r+1,j+1)*Xkj{r};
            end
            
            
            
            switch js
                case 1
                    Upar = vertcat(Uk,p_case(1));
                case 2
                    if k ==0
                        Upar = vertcat(Uk,p_case(1));
                    else
                        Upar = vertcat(Uk,p_case(2));
                    end
                case 3
                    if k ==0
                        Upar = vertcat(Uk,p_case(1));
                    else
                        Upar = vertcat(Uk,p_case(3));
                    end
                case 4
                    if k ==0
                        Upar = vertcat(Uk,p_case(2));
                    else
                        Upar = vertcat(Uk,p_case(1));
                    end
                case 5
                    Upar = vertcat(Uk,p_case(2));
                case 6
                    if k ==0
                        Upar = vertcat(Uk,p_case(2));
                    else
                        Upar = vertcat(Uk,p_case(3));
                    end
                case 7
                    if k ==0
                        Upar = vertcat(Uk,p_case(3));
                    else
                        Upar = vertcat(Uk,p_case(1));
                    end
                case 8
                    if k ==0
                        Upar = vertcat(Uk,p_case(3));
                    else
                        Upar = vertcat(Uk,p_case(2));
                    end
                case 9
                    Upar = vertcat(Uk,p_case(3));
            end
            
            
            [fj,qj] =  f(Xkj{j},Upar);
            
            g = {g{:},par.tf*fj-xp};  % dynamics and algebraic constraints
            lbg = [lbg;zeros(nx,1)];
            ubg = [ubg;zeros(nx,1)];
            
            % Add contribution to the end states
            Xk_end = Xk_end + D(j+1)*Xkj{j};
            
            J = J + (B(j+1)*qj*par.tf );
            
        end
        
        
        Uk_prev = MX.sym(['Uprev_' num2str(k+1) '_' num2str(js)],nu);
        Uk_prev = Uk;
        
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
        
        
        if k < nR
            u_ant{js,k+1} = MX.sym(['u_ant_' num2str(k) '_' num2str(js)],nu);
            u_ant{js,k+1} = Uk;
        end
        
    end
end

%Add Non-anticipativity constraints
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

opts = struct('warn_initial_bounds',false, ...
    'print_time',false, ...
    'ipopt',struct('print_level',1) ...
    );

nlp = struct('x',vertcat(w{:}),'p',vertcat(p,U0,X0_par),'f',J,'g',vertcat(g{:}));

solver = nlpsol('solver','ipopt',nlp,opts);

%% simulation
u_in = u0;
for sim_k = 1:par.nIter
    
    nlp_par = vertcat(p_case,u_in,dx0);
    tic
    sol = solver('x0',w0,'p',nlp_par,'lbx',lbw,'ubx',ubw,'lbg',lbg,'ubg',ubg);
    sim.t_solve(sim_k) = toc;
    
    flag = solver.stats();
    disp(['Iter ' num2str(sim_k) ': ' flag.return_status])
    
    % extract solution
    w_opt = full(sol.x);
    n_w_i = nx+par.N*(nu+(d+1)*nx);
    for i = 1:nS
        u_opt(:,i) = [w_opt((i-1)*n_w_i+nx+1:nu+(d+1)*nx:i*n_w_i);NaN];
        x1_opt(:,i) = w_opt([(i-1)*n_w_i+1,(i-1)*n_w_i+nx+nu+1:nu+(d+1)*nx:i*n_w_i]);
        x2_opt(:,i) = w_opt([(i-1)*n_w_i+2,(i-1)*n_w_i+nx+nu+2:nu+(d+1)*nx:i*n_w_i]);
    end
    
    % inject first control input into the plant
    u_in = u_opt(1,1);
    
    % plant simulator
    Fk = F('x0',dx0,'p',[u_in;p_val]);
    dx0 = full(Fk.xf);
    
    sim.u(sim_k) = u_in;
    sim.x(sim_k,:) = dx0;
    
    % plot closed loop implementation and open loop predictions
    figure(12)
    clf
    subplot(211)
    hold all
    plot(-sim_k+1:0,sim.x(:,2))
    plot(0:par.N,x2_opt)
    xlim([-par.nIter,par.N])
    subplot(212)
    hold all
    stairs(-sim_k+1:0,sim.u)
    stairs(0:par.N,u_opt)
    xlim([-par.nIter,par.N])
    
    % re-initialize upper and lower bounds
    lbg = [];
    ubg = [];
    for js = 1:nS
        lbg = [lbg;dx0];
        ubg = [ubg;dx0];
        for k = 0:par.N-1
            for j = 1:d
                lbg = [lbg;zeros(nx,1)];
                ubg = [ubg;zeros(nx,1)];
            end
            lbg = [lbg;zeros(nx,1)];
            ubg = [ubg;zeros(nx,1)];
        end
    end
    for k = 1:nR
        for js = 1:nS-1
            lbg = [lbg;zeros(nu,1)];
            ubg = [ubg;zeros(nu,1)];
        end
    end
    
    
end

figure(11)
subplot(211)
hold all
plot(sim.x)
subplot(212)
hold all
stairs(sim.u)

%% save result

choice = questdlg('Save results?', ...
    'Save Result', ...
    'Yes','No','No');
% Handle response
switch choice
    case 'Yes'
        save('sim_central3.mat','sim')
end