function [t,states,xdot,inputs] = cstrFunc()
%CSTRFUNC Summary of this function goes here
% 
% [OUTPUTARGS] = CSTRFUNC(INPUTARGS) Explain usage here
% 
% Examples: 
% 
% Provide sample usage code here
% 
% See also: List related files here

% $Author: suwartad $	$Date: 2017/10/07 16:38:51 $	$Revision: 0.1 $
% Copyright: Process Control Group - NTNU Trondheim 2017

import casadi.* 

V   = 10;   % Volume in [L]
k   = 1.2;  % Reaction rate in [L / (mol * minute)]
Caf = 1;    % Feeding concentration [mol/L]

% symbolic primitives
t  = SX.sym('t');
x1 = SX.sym('x1');  % Ca
x2 = SX.sym('x2');  % Cb
states = [x1;x2];
Q      = SX.sym('Q');
inputs = Q;

% CSTR model
x1dot  = (Q/V)*(Caf - x1) - k*x1;
x2dot  = (Q/V)*-x2 + k*x1;
xdot   = [x1dot;x2dot];

end
