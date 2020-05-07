%% Implementation of a Gauss-Neston SQP solver usign CasADi
close all
clear variables
clc

% import CasADi
import casadi.*

nv = 2;
x = MX.sym('x',nv);

x_test = [1.1 3.2];
% YOUR CODE HERE: define the objective function (Task 5.2)
% f     = ...;
% F     = ...;
% Jf    = ...;
f  = 0.5*(x(1) - 1)^2 + 0.5*(10*(x(2) - x(1)^2))^2 + 0.5*x(2)^2;
F  = Function('F',{x},{f});
Jf = Function('Jf',{x},{jacobian(f,x)});
F_test  = 203.1300;
Jf_test = [ -437.7000,  202.2000 ];

% YOUR CODE HERE: define the residuals (Task 5.2)
% r     = ...;
% R     = ...;
% Jr    = ...;
r  = [x(1)-1;10*(x(2)-x(1)^2);x(2)];
R  = Function('R',{x},{r});
Jr = Function('Jr',{x},{jacobian(r,x)});
R_test = [ 0.1000; 19.9000; 3.2000 ];
Jr_test = [ 1,     0; -22,    10;  0,     1 ];

% YOUR CODE HERE: define the equality constraints 
% g     = ...;
% G     = ...;
% Jg    = ...;
g  = x(1) + (1 - x(2))^2;
G  = Function('G',{x},{g});
Jg = Function('Jg',{x},{jacobian(g,x)});
G_test  =  5.9400;
Jg_test = [ 1.0000,    4.4000];

% YOUR CODE HERE: define the inequality constraints 
% h     = ...;
% H     = ...;
% Jh    = ...;
h   = 0.2 + x(1)^2 - x(2);
H   = Function('H',{x},{h});
Jh  = Function('Jh',{x},{jacobian(h,x)});
H_test =  -1.7900;
Jh_test =  [ 2.2000,   -1.0000 ];

% Define linearization point
xk = MX.sym('xk',nv);

% YOUR CODE HERE: define the linearized equalities 
% Jtemp     = Jg({xk});
% g_temp    = G({xk});
% g_l       = ...;
Jgtemp    = Jg(xk);
g_temp    = G(xk);
g_l       = g_temp + Jgtemp*x;

% YOUR CODE HERE: define the linearized inequalities 
% Jtemp     = Jh({xk});
% g_temp    = H({xk});
% h_l       = ...;
Jhtemp    = Jh(xk);
h_temp    = H(xk);
h_l       = h_temp + Jhtemp*x;

% YOUR CODE HERE: Gauss-Newton Hessian approximation (Task 5.3)
% j_out     = Jr({xk});
% jf_out    = Jf({xk});
% r_out     = R({xk});
% f_gn      = ...;
j_out     = Jr(xk);
jf_out    = Jf(xk);
r_out     = R(xk);
f_gn      = 0.5*(x'*(j_out)'*j_out*x) + jf_out*x;

% Allocate QP solver
% qp = struct('x',x, 'f',f_gn,'p',xk);
% solver = qpsol('solver', 'qpoases', qp);

% qp = struct('x',x, 'f',f_gn,'g',g_l,'p',xk);
% solver = qpsol('solver', 'qpoases', qp);

qp = struct('x',x, 'f',f_gn,'g',[g_l;h_l],'p',xk);
solver = qpsol('solver', 'qpoases', qp);

% SQP solver
max_it = 100;
xk = [1 1]'; % Initial guess
iter = zeros(nv,max_it);
iter(:,1) = xk;

for i=2:max_it    
    % YOUR CODE HERE: formulate the QP (Tasks 5.3, 5.4, 5.5)
    % arg     = struct;
    % arg.lbg =  ...;           
    % arg.ubg =  ...;
    % arg.lbx =  ...;           
    % arg.ubx =  ...; 
    % arg.p   =  ...;
    
%     arg     = struct;
%     %arg.lbg =  0;           
%     %arg.ubg =  0;
%     arg.lbg =  [0;-inf];           
%     arg.ubg =  [0;0];
%     arg.lbx =  -inf;           
%     arg.ubx =  inf; 
%     arg.p   =  xk;

    % Solve with qpOASES
    %sol = solver(arg);
    sol  = solver('lbg',[0;-inf],'ubg',[0;0],'lbx',-inf,'ubx',inf,'p',xk);
    step = full(sol.x);
    if norm(step) < 1.0e-16
        break;
    end
 
    t = 1;
    
    iter(:,i) = iter(:,i-1) + t*step;
    xk = iter(:,i);
end

% Plot results
iter = iter(:,1:i-1);
[X,Y] = meshgrid(-1.5:.05:1.5, -1.5:.05:1.5);
Z = log(1 + 1/2*(X -1).^2 + 1/2*(10*(Y -X.^2)).^2 + 1/2*Y.^2);
figure()
subplot(1,2,1)
surf(X,Y,Z)
xlabel('x_1')
ylabel('x_2')
hold all
plot3(iter(1,:),iter(2,:),zeros(length(iter(1,:)),1),'black')
plot3(iter(1,:),iter(2,:),zeros(length(iter(1,:)),1),'blacko')
y_g = linspace(-0.08,1.1,20);
x_g = -(1 - y_g).^2;
plot3(x_g,y_g,zeros(length(x_g),1),'r');

x_h = linspace(-1.5,1.5,20);
y_h = 0.2 + x_h.^2;
plot3(x_h,y_h,zeros(length(x_h),1),'r--');

xlim([-1.5,1.5]);
ylim([-1.5,1.5]);
subplot(1,2,2)
contour(X,Y,Z)
hold all
plot3(iter(1,:),iter(2,:),zeros(length(iter(1,:)),1),'black')
plot3(iter(1,:),iter(2,:),zeros(length(iter(1,:)),1),'blacko')
xlabel('x_1')
ylabel('x_2')
plot3(x_g,y_g,zeros(length(x_g),1),'r');
plot3(x_h,y_h,zeros(length(x_h),1),'r--');
xlim([-1.5,1.5]);
ylim([-1.5,1.5]);
figure()
plot(iter(1:2,:)')
grid on
xlabel('iterations')
ylabel('primal solution')
grid on
