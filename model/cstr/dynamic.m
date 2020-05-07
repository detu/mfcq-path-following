function [x, t_intermediate, x_intermediate] = dynamic(system, T, t0, x0, u0)
%DYNAMIC Summary of this function goes here
% 
% [OUTPUTARGS] = DYNAMIC(INPUTARGS) Explain usage here
% 
% Examples: 
% 
% Provide sample usage code here
% 
% See also: List related files here

% $Author: suwartad $	$Date: 2018/04/19 05:00:34 $	$Revision: 0.1 $
% Copyright: Process Control Group - NTNU Trondheim 2018

x = system(t0, x0, u0, T);
x_intermediate = [x0; x];
t_intermediate = [t0, t0+T];

end
