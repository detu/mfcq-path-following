function [Eta, z] = computeEta(Jeq, g, y, cin, activeBoundTol)
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
%z       = g - Jeq'*y.lam_g;
%z       = g + Jeq'*y.lam_g;

% active bound constraints
boundMultiplier    = y.lam_x;
positiveBoundMult  = abs(boundMultiplier);
if(isempty(activeBoundTol))
    activeBoundTol = 1e-1;
end
activeIndex        = find(positiveBoundMult>1e-1);  % set active bound constraint.
% activeIndex        = find(positiveBoundMult>activeBoundTol);
paramIndex         = find(activeIndex <= 84);       % remove the first 84 constraints (the parameter)
activeIndex(paramIndex) = [];
numActiveBoundCons = numel(activeIndex);

[numY,m] = size(Jeq);
JAbc     = zeros(numActiveBoundCons,m);
for i=1:numActiveBoundCons
    row         = activeIndex(i);
    JAbc(i,row) = 1;
end

Jeq     = [Jeq;JAbc];

%z       = g + Jeq'*y.lam_g + y.lam_x;
z       = g + Jeq'*([y.lam_g; y.lam_x(activeIndex)]);
stackLC = [cin;z];
Eta     = 0.1*norm(stackLC,inf);

end
