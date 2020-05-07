function [u_nlp_opt,plotState] = plotStatesSoftConstraint(data, N)
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
nsc = 2;

% optimized initial state
x0_opt       = data(1:nx);  
data(1:nx)   = [];
data         = reshape(data, (nu + (nx+nsc+ns)*d + (nx+nsc+ns)), N*nk);
u_nlp_opt    = data(1:nu,1:N*nk);
data(1:nu,:) = [];     % remove optimized controls

% lb0          = lb(1:nx+ns);
% lb(1:nx)     = [];
% lb           = reshape(lb, (nu + (nx+nsc+ns)*d + (nx+nsc+ns)), N*nk);
% lbU          = lb(1:nu,1:N*nk);
% lb(1:nu,:)   = [];
% 
% ub0          = ub(1:nx+ns);
% ub(1:nx)     = [];
% ub           = reshape(ub, (nu + (nx+nsc+ns)*d + (nx+nsc+ns)), N*nk);
% ubU          = ub(1:nu,1:N*nk);
% ub(1:nu,:)   = [];


% prepare a matrix for plotting
%nState    = (nx+ns) + N*nk*(d+1)*(nx+ns);
nState    = (nx+nsc+ns) + N*nk*(d+1)*(nx+nsc+ns);
nPoint    = nState/(nx+nsc+ns);
plotState = zeros(nx+nsc+ns,nPoint);
plotState(1:nx,1) = x0_opt;

% plotLb         = zeros(nx+ns,nPoint);
% plotLb(:,1)    = lb0;
% 
% plotUb         = zeros(nx+ns,nPoint);
% plotUb(:,1)    = ub0;

% extract states from each collocation point and each time horizon
sInd    = 2; % initial index row
for i=1:N*nk
    temp    = data(:,i);
    numCol  = size(temp,1);
    numRow  = numCol/(nx+nsc+ns);
    temp    = reshape(temp,nx+nsc+ns,numRow);
    plotState(:,sInd:(numRow+sInd-1)) = temp;
    
%     tempLb = lb(:,i);
%     tempLb = reshape(tempLb,nx+nsc+ns,numRow);
%     plotLb(:,sInd:(numRow+sInd-1)) = tempLb;
%     
%     tempUb = ub(:,i);
%     tempUb = reshape(tempUb,nx+nsc+ns,numRow);
%     plotUb(:,sInd:(numRow+sInd-1)) = tempUb;
    
    sInd    = numRow + sInd;
end

end
