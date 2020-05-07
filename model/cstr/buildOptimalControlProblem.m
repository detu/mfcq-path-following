function [J,g,w0,w,lbg,ubg,lbw,ubw,params] = buildOptimalControlProblem(optProblem, N, x0, u0, x0_measure)
%BUILDOPTIMALCONTROLPROBLEM Summary of this function goes here
% 
% [OUTPUTARGS] = BUILDOPTIMALCONTROLPROBLEM(INPUTARGS) Explain usage here
% 
% Examples: 
% 
% Provide sample usage code here
% 
% See also: List related files here

% $Author: suwartad $	$Date: 2018/04/19 05:31:40 $	$Revision: 0.1 $
% Copyright: Process Control Group - NTNU Trondheim 2018
import casadi.*
global nx;
%% supply initial guess
% call ODE15s N-times initial guess in optimization
x(1,:) = x0(1:nx)';
for k=1:N
    x(k+1,:) = x0(1:nx)';    % ONE-TIME SIMULATION !
end

%% construct variables for a realization
[J,g,w0,w,lbg,ubg,lbw,ubw,params] = optProblem(x, u0, N, x0_measure);



end
