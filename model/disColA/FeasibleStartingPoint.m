%FEASIBLESTARTINGPOINT Summary of this function goes here
%
% Finding a feasible starting point for collocation points
%
% [OUTPUTARGS] = FEASIBLESTARTINGPOINT(INPUTARGS) Explain usage here
% 
% Examples: 
% 
% Provide sample usage code here
% 
% See also: List related files here

% $Author: suwartad $	$Date: 2019/02/08 14:26:47 $	$Revision: 0.1 $
% Copyright: Process Control Group - NTNU Trondheim 2019

import casadi.*

numberOfCollocationPoint  = 3;
%numberOfPredictiveHorizon = 30;
numberOfPredictiveHorizon = 90;
nx = 84; % number of state variable
nu = 5;  % number of control variable


%% collocation routine
% Degree of interpolating polynomial
d = numberOfCollocationPoint;

% Get collocation points
tau_root = [0 collocation_points(d, 'legendre')];

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

%% model part

% Time horizon
T = numberOfPredictiveHorizon;

NT = 41;
Uf = 0.31;           % Feeding rate F_0
%Uf = 0.30;

% invoke the model
[~,state,xdot,inputs] = DistColACstr(Uf);

% objective function value
obj = computeObjectiveFunction(state,inputs);

f = Function('f',{state,inputs}, {xdot,obj});

%% optimization part
% Control discretization
N = numberOfPredictiveHorizon; % number of control intervals
h = T/N;

% Start with an empty NLP
w   = {};
w0  = [];
lbw = [];
ubw = [];
J   = 0;
g   = {};
lbg = [];
ubg = [];

% load starting point with Feeding Rate 0.31
%load Xinit31.mat;
%x0 = Xinit31(1:nx);
%u0 = Xinit31(85:89);
% load Xinit30.mat
% x0 = Xinit30(1:nx);
%u0 = Xinit30(85:89);
% load Xinit3075.mat
% x0 = Xinit3075(1:nx);
%u0 = Xinit3075(85:89);

%load Xinit30Ub65.mat
%x0 = Xinit30Ub65(1:nx);
%u0 = Xinit30Ub65(85:89);
load Xinit31Ub64.mat
x0 = Xinit31Ub64(1:nx);
u0 = Xinit31Ub64(85:89);


% bound constraints
VB_max = 4.008;
xB_max = 0.1;   

% State bounds and initial guess
x_min =  zeros(84,1);  % try without epsilon here, later put epsilon
x_max =  ones(84,1);

% x_min(84) = 0.3;
% x_max(1)  = xB_max; % scaled bottom concentration
x_max(84) = 0.65;
%x_max(84) = 0.75;
x_min(43) = 0.5;
x_max(43) = 0.5;
x_min(83) = 0.5;
x_max(83) = 0.5;

% Set fixed control bounds
% load Xinit31.mat;
% u0 = Xinit31(85:89);
% load Xinit3175.mat;
% u0 = Xinit3175(85:89);
% u5 = Xinit31(89);

% load Xinit31Ub65.mat;
% u0 = Xinit31Ub65(85:89);
u_min = [u0(1); u0(2); u0(3); u0(4); u0(5)];
u_max = [u0(1); u0(2); u0(3); u0(4); u0(5)];
% u_min = [0; 0; 0; 0; 0];
% u_max = [inf; inf; inf; inf; inf];

Xopt = [];

% Formulate the NLP
for k=0:N-1
    
    % "Lift" initial conditions -x0 should be update from previous end
    % point of collocation
    Xk  = MX.sym('X0', nx);
    w   = {w{:}, Xk};
    lbw = [lbw; x_min];
    ubw = [ubw; x_max];
    w0  = [w0; x0];
    g   = {g{:}, Xk - x0};
    lbg = [lbg; zeros(nx,1)];
    ubg = [ubg; zeros(nx,1)];
    
    
    % New NLP variable for the control
    Uk     = MX.sym(['U_' num2str(k)], nu);
    w      = {w{:}, Uk};
    lbw    = [lbw; u_min];
    ubw    = [ubw; u_max];
    w0     = [w0;  u0];

    % State at collocation points
    Xkj = {};
    for j=1:d
        Xkj{j} = MX.sym(['X_' num2str(k) '_' num2str(j)], nx);
        w      = {w{:}, Xkj{j}};
        lbw    = [lbw; x_min];
        ubw    = [ubw; x_max];
        w0     = [w0; x0];
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
       [fj, qj] = f(Xkj{j},Uk);
       g   = {g{:}, h*fj - xp};
       lbg = [lbg; zeros(nx,1)];
       ubg = [ubg; zeros(nx,1)];

       % Add contribution to the end state
       Xk_end = Xk_end + D(j+1)*Xkj{j};

       % Add contribution to quadrature function
       %J = J + B(j+1)*qj*h;
       J = 1;
    end

    % New NLP variable for state at end of interval
    Xk  = MX.sym(['X_' num2str(k+1)], nx);
    w   = {w{:}, Xk};
    lbw = [lbw; x_min];
    ubw = [ubw; x_max];
    w0  = [w0; x0];
    
    % Add equality constraint
    g   = {g{:}, (Xk_end-Xk)};
    lbg = [lbg; zeros(nx,1)];
    ubg = [ubg; zeros(nx,1)];

    % call optimizer
    prob   = struct('f', J, 'x', vertcat(w{:}), 'g', vertcat(g{:}));
    solver = nlpsol('solver', 'ipopt', prob);
    
    % Solve the NLP
    sol   = solver('x0', w0, 'lbx', lbw, 'ubx', ubw,'lbg', lbg, 'ubg', ubg);
    w_opt = full(sol.x);
    
    success = strcmp(solver.stats.return_status,'Infeasible_Problem_Detected');
    if (success)
        keyboard;
    end
    
    % arrange result
    XkOpt = ArrangeFeasibleOptResults(w_opt,d);
    
    % copy the last optimized Xk_end
    x0 = XkOpt(:,d+1);
    
    % append the whole XkOpt for all collocation point
    Xopt = [Xopt XkOpt];

end

keyboard;