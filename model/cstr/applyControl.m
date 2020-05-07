function [tapplied, xapplied] = applyControl(system, T, t0, x0, u0)
%APPLYCONTROL Summary of this function goes here
% 
% [OUTPUTARGS] = APPLYCONTROL(INPUTARGS) Explain usage here
% 
% Examples: 
% 
% Provide sample usage code here
% 
% See also: List related files here

% $Author: suwartad $	$Date: 2018/04/19 04:53:26 $	$Revision: 0.1 $
% Copyright: Process Control Group - NTNU Trondheim 2018

xapplied = dynamic(system, T, t0, x0, u0(:,1));
tapplied = t0+T;

end
