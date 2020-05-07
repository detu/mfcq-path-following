function xOpt = ArrangeFeasibleOptResults(w,d)
%ARRANGEFEASIBLEOPTRESULTS Summary of this function goes here
% 
% [OUTPUTARGS] = ARRANGEFEASIBLEOPTRESULTS(INPUTARGS) Explain usage here
% 
% Examples: 
% 
% Provide sample usage code here
% 
% See also: List related files here

% $Author: suwartad $	$Date: 2019/02/08 16:16:52 $	$Revision: 0.1 $
% Copyright: Process Control Group - NTNU Trondheim 2019

nx = 84;
nu = 5;

% remove initial state and control
w(1:nx+nu) = [];

xOpt  = zeros(nx,d+1);
indexStart = 1;
indexEnd   = nx;
for i=1:d+1
    xOpt(:,i)   = w(indexStart:indexEnd);
    indexStart  = indexEnd + 1;
    indexEnd    = indexStart + nx - 1;
end


end
