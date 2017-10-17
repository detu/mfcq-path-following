function cstrSsOpt
%CSTRSSOPT Summary of this function goes here
% 
% Steady-state optimization of an isothermal CSTR
% reaction: A -> B
% state variables (controlled variable): Ca and Cb
% control input (manipulated variable): Q
%
% [OUTPUTARGS] = CSTRSSOPT(INPUTARGS) Explain usage here
% 
% Examples: 
% 
% Provide sample usage code here
% 
% See also: List related files here

% $Author: suwartad $	$Date: 2017/10/06 22:29:36 $	$Revision: 0.1 $
% Copyright: Process Control Group - NTNU Trondheim 2017

format long;
import casadi.*

% parameter values
V   = 10;   % Volume in [L]
k   = 1.2;  % Reaction rate in [L / (mol * minute)]
Caf = 1;    % Feeding concentration [mol/L]

% symbolic primitives
u  = SX.sym('u');   % Q
x1 = SX.sym('x1');  % Ca
x2 = SX.sym('x2');  % Cb

% concatenate states and controls 
x   = [x1;x2;u];

% initial guess
Uinit = [0.1;0.1;15];

% define the dynamics as equality constraints and additional inequality constraints
[obj,eq, lbx, ubx, lbg, ubg] = buildModelEq(x,V,k,Caf);

prob    = struct('f', obj, 'x', x, 'g', eq);
options = struct;
%options.ipopt.tol       = 1e-12;
%options.acceptable_compl_inf_tol    = 1e-6;
solver  = nlpsol('solver', 'ipopt', prob, options);

% Solve the NLP
startnlp = tic;
sol   = solver('x0', Uinit, 'lbx', lbx, 'ubx', ubx, 'lbg', lbg, 'ubg', ubg);
elapsednlp = toc(startnlp);
fprintf('IPOPT solver runtime = %f\n',elapsednlp);

u      = full(sol.x);
lamda  = full(sol.lam_g);
Xinit  = u;
save CstrXinit.mat Xinit;
save LamdaCstr.mat lamda;

%% Compute Hessian and perform Greshgorin convexification

% symbolic variable for dual variables
l1 = SX.sym('l1');  
l2 = SX.sym('l2');  

xsol = u;
lambda.eqnonlin = lamda;
% extract Lagrange multiplier
%l   = vertcat(l{:});
l  = [l1;l2];
L  = obj + l'*eq;

Lagr = Function('Lagr', {x,l}, {L}, char('x','l'), char('Lagr'));
%Jac  = Function(Lagr.jacobian('x','Lagr'));
H    = Function(Lagr.hessian('x','Lagr'));
cons = Function('Const', {x}, {eq}, char('x'), char('cons'));
Jcon = Function(cons.jacobian('x','cons'));

eqVal = cons(xsol);
eqVal = full(eqVal);
Hx   = H(xsol,lambda.eqnonlin);
Hx   = full(Hx);
Jac  = Jcon(xsol);
Jac  = full(Jac);
rH   = null(Jac)'*Hx*null(Jac)
erH  = eig(rH)


[Hxxl,Qmax]   = Greshgorin(Hx);
save QmaxCstr.mat Qmax;

% check at initial point for optimization
xstat = Xinit(1:2);
u0    = 19;
xeval = [xstat;u0];
Jeval  = Jcon(xeval);
Jeval  = full(Jeval);
Hxxl   = H(xeval,lambda.eqnonlin);
Hxxl   = full(Hxxl);
Hconv  = Hxxl + diag(Qmax);
rHe    = null(Jeval)'*Hconv*null(Jeval);


keyboard;
end

function [J, ceq, lbx, ubx, lbg, ubg] = buildModelEq(u,V,k,Caf)
import casadi.* 

Ca = u(1);
Cb = u(2);
Q  = u(3);

% objective function 
J = -Q*(2*Cb - 0.5);

% CSTR model
% define xdot
x1dot  = SX.sym('x1dot');
x2dot  = SX.sym('x1dot');

x1dot  = (Q/V)*(Caf - Ca) - k*Ca;
x2dot  = (Q/V)*-Cb + k*Ca;
ceq    = [x1dot;x2dot];

% bound constraints
lbx  = [0;0;10];
ubx  = [1;1;20];
lbg  = [0;0];
ubg  = [0;0];

end

function [H,Q] = Greshgorin(H)
numH    = size(H,1);
Q       = zeros(numH,numH);
%Q       = eye(numH);
%delta   = 1e-1;
%delta   = 1e-2;
%delta   = 2.5e-1;  % normal case
%delta   = 0.5;
%delta   = 1;
%delta   = 1.5;
%delta   = 2;
%delta   = 2.5;      % with measurement noise 1 percent
%delta   = 5;
%delta   = 100;
delta   = 1e-6;
for i=1:numH  % iterate all row of Hessian
    sumRow = 0;
    for j=1:numH
        if j ~= i
            sumRow = sumRow + abs(H(i,j));
        end
    end
    %sumRow = sumRow - abs(H(i,i));
    
    if H(i,i) <= sumRow   % include equality 
    %if H(i,i) < sumRow 
        %Q(i,i) = sumRow - H(i,i) + delta;
        Q(i,i) = sumRow - H(i,i) + delta;
    end
end

% % loop over Qm to obtain maximum number
% for i=1:numH
%     if Q(i,i) > Qm(i,1)
%         Qm(i,1) = Q(i,i);
%     end
% end
Q = diag(Q);
end

