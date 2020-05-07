function [u_nlp_opt,plotState] = plotStatesDirectMS(data, N)
%PLOTSTATESN Summary of this function goes here
% 
% [OUTPUTARGS] = PLOTSTATESN(INPUTARGS) Explain usage here
% 
% Examples: 
% 
% Provide sample usage code here
% 
% See also: List related files here

% $Author: suwartad $	$Date: 2016/04/27 22:42:29 $	$Revision: 0.1 $
% Copyright: Process Control Group - NTNU Trondheim 2016

global nk nx nu d ns;

% include soft constraint variables to the states
nsc = 0;
d = 0;

% optimized initial state
x0_opt       = data(1:nx);  
data(1:nx)   = [];
data         = reshape(data, (nu + (nx+nsc+ns)*d + (nx+nsc+ns)), N*nk);
u_nlp_opt    = data(1:nu,1:N*nk);
data(1:nu,:) = [];     % remove optimized controls


% prepare a matrix for plotting
nState    = (nx+nsc+ns) + N*nk*(d+1)*(nx+nsc+ns);
nPoint    = nState/(nx+nsc+ns);
plotState = zeros(nx+nsc+ns,nPoint);
plotState(1:nx,1) = x0_opt;


% extract states from each collocation point and each time horizon
sInd    = 2; % initial index row
for i=1:N*nk
    temp    = data(:,i);
    numCol  = size(temp,1);
    numRow  = numCol/(nx+nsc+ns);
    temp    = reshape(temp,nx+nsc+ns,numRow);
    plotState(:,sInd:(numRow+sInd-1)) = temp;
    
    sInd    = numRow + sInd;
end

end
