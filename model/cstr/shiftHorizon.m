function u0 = shiftHorizon(u)
%SHIFTHORIZON Summary of this function goes here
% 
% [OUTPUTARGS] = SHIFTHORIZON(INPUTARGS) Explain usage here
% 
% Examples: 
% 
% Provide sample usage code here
% 
% See also: List related files here

% $Author: suwartad $	$Date: 2018/04/19 04:55:54 $	$Revision: 0.1 $
% Copyright: Process Control Group - NTNU Trondheim 2018

u0 = [u(:,2:size(u,2)) u(:,size(u,2))];

end
