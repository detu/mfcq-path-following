function [Eta, z] = computeEta(Jeq, g, y, cin)
%COMPUTEETA Summary of this function goes here
% 
% [OUTPUTARGS] = COMPUTEETA(INPUTARGS) Explain usage here
% 
% Examples: 
% 
% Provide sample usage code here
% 
% See also: List related files here

% $Author: suwartad $	$Date: 2017/05/18 20:20:52 $	$Revision: 0.1 $
% Copyright: Process Control Group - NTNU Trondheim 2017

% only equality constraints 
% compute Lagrangian
z       = g - Jeq'*y.lam_g;
stackLC = [cin;z];
Eta     = norm(stackLC,inf);

end
