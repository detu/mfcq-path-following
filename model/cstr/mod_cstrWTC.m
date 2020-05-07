function xprime=mod_cstrWTC(t,X,U) 
%MOD_CSTR Summary of this function goes here
% 
% [OUTPUTARGS] = MOD_CSTR(INPUTARGS) Explain usage here
% 
% Examples: 
% 
% Provide sample usage code here
% 
% See also: List related files here

% $Author: suwartad $	$Date: 2017/10/07 16:33:12 $	$Revision: 0.1 $
% Copyright: Process Control Group - NTNU Trondheim 2017


%% CSTR model
%global k;
% parameter values
V   = 10;   % Volume in [L]
k   = 1.2;  % Reaction rate in [L / (mol * minute)]
Caf = 1;    % Feeding concentration [mol/L]

% state equations
Q      = U;  
Ca     = X(1);
Cb     = X(2);
x1dot  = (Q/V)*(Caf - Ca) - k*Ca;
x2dot  = (Q/V)*-Cb + k*Ca;

% Output
xprime=[x1dot;x2dot];
end
