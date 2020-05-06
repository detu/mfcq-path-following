function L = computeObjectiveFunction(x,u)
%COMPUTEOBJECTIVEFUNCTION Summary of this function goes here
% 
% [OUTPUTARGS] = COMPUTEOBJECTIVEFUNCTION(INPUTARGS) Explain usage here
% 
% Examples: 
% 
% Provide sample usage code here
% 
% See also: List related files here

% $Author: suwartad $	$Date: 2019/02/06 13:29:59 $	$Revision: 0.1 $
% Copyright: Process Control Group - NTNU Trondheim 2019

global nx nu;
nx = 84;
nu = 5;
load CstrDistXinit.mat;
xf    = Xinit(1:84);
u_opt = Xinit(85:89);

load Qmax.mat;

F_0 = 0.3;
% prices
pf = 1;
pV = 0.02;
pB = 2;
pD = 0;


Jstate     = (Qmax(1:nx,1).*(x - xf))' * (x - xf);
Jecon      = (pf*F_0 + pV*u(2) - pB*u(5));
L          = Jecon + Jstate;

end
