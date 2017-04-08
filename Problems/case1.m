function [prob] = case1(p)
% from V.K and J.J paper
% PathfollowDG(@(p)case1(p), 0,1,[0;0;0],[1e-8;0.3333;0.3333;0.3333;4e-9;1.6e-9;2.6e-9],1.2)
% PathfollowDG(@(p)case1(p), 0,1,[0;0;0],[-7.52941170969587e-09;0.3333;0.3333;0.3333;3.999999999999999e-09;1.599999992320000e-09;2.666666645333332e-09],1.2)
% PathfollowDG(@(p)case1(p), 0,1,[0;0;0],[-7.52941170969587e-09;0;1;0;3.999999999999999e-09;1.599999992320000e-09;2.666666645333332e-09],1.2)
% PathfollowDG(@(p)case1(p), 0,1,[0;0;0],[0;0.5;0;0.5;0;0;0],1.2)
% PathfollowDG(@(p)case1(p), 0,1,[0;0;0],[0;0;1;0;0;0;0],1.2)
prob = struct('n',{0},'m',{0},'bl',{0},'bu',{0},'cl',{0},'cu',{0},'obj',{0},'cons',{0},'hess',{0},'x',0,'name',0);
prob.n = 3;
prob.m = 7;
prob.me = 1;
%prob.cl = [0;p(2);p(1)];
%prob.cu = [Inf;Inf;p(1)];

prob.dcdp = [-10; 0; 10; 20; 0; -10; -10];

prob.name = 'Problem 1';
prob.x = zeros(3,1);

%prob.obj  = (@(x,p)(objective(x,p)));
prob.obj  = (@(x)(objective(x,p)));
%prob.cons = (@(x,p)(constraint(x,p)));
prob.cons = (@(x)(constraint(x,p)));
%prob.hess = (@(x,y,p)(hessian(x,y,p)));
prob.hess = (@(x,y)(hessian(x,y,p)));
%prob.hc   = (@(x,y,p)(hessianConstraints(x,y,p)));
prob.chess = (@(x)(conhess(x,p)));
end

function [f,g] = objective(x,p)

f = -exp(x(2)) + 0.5*(x(1) - x(3))^2;
g = zeros(3,1);
g(1) =  x(1) - x(3);
g(2) = -exp(x(2));
g(3) = -x(1) + x(3);
end

function [c,J] = constraint(x,p)
c = [x(3) - 10*p; ...
     x(1) - x(2); ...
     10*p - x(2); ...
     -x(1) - x(2) + 20*p; ...
     5 - x(1); ...
     0.5*x(1) - x(2) + 7.5 - 10*p; ...
%     0.5*x(1) - x(2) + 12.5 - 10*p; ...
%     -0.5*x(1) - x(2) + 7.5 - 10*p];
    -0.5*x(1) - x(2) + 12.5 - 10*p];

 J = zeros(7,3);
 J(1,3) = 1;
 J(2,1) = 1;
 J(2,2) = -1;
 J(3,2) = -1;
 J(4,1) = -1;
 J(4,2) = -1;
 J(5,1) = -1;
 J(6,1) = 0.5;
 J(6,2) = -1;
 J(7,1) = -0.5;
 J(7,2) = -1;
end

function [H] = hessian(x,y,p)
H = zeros(3,3);
y = y;
H(1,1) = 1;
H(1,3) = -1;
H(2,2) = -exp(x(2));
H(3,1) = -1;
H(3,3) = 1;
end

% function [Hc] = hessianConstraints(x,y,p)
% Hc = zeros(7,7);
% end

function HC = conhess(x,p)
HC =  zeros(3,3,7);
end

