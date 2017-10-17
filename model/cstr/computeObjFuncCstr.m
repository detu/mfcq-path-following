function Jobj = computeObjFuncCstr(uOpt,xActual)
%COMPUTEOBJFUNCCSTR Summary of this function goes here
% 
% [OUTPUTARGS] = COMPUTEOBJFUNCCSTR(INPUTARGS) Explain usage here
% 
% Examples: 
% 
% Provide sample usage code here
% 
% See also: List related files here

% $Author: suwartad $	$Date: 2017/10/07 21:52:18 $	$Revision: 0.1 $
% Copyright: Process Control Group - NTNU Trondheim 2017


% steady-state values
load CstrXinit.mat;
xs    = Xinit(1:2);
us    = Xinit(3);
nx    = size(xs,1);
nu    = size(us,1);
% weights
load QmaxCstr.mat;
c1  = -0.1; % noisy case
lss = -5.999999999806272 + c1; %steady-state objective function value

Jecon    = -uOpt*(2*xActual(2) - 0.5);
Jcontrol = (Qmax(nx+1:nx+nu,1).*(uOpt - us))' * (uOpt - us);
Jstate   = (Qmax(1:nx,1).*(xActual - xs))' * (xActual - xs);%
J        = Jecon + Jcontrol + Jstate - lss;
%J        = Jecon + Jcontrol + Jstate;
fprintf('-------------------------------------------------------------------- \n');
fprintf('Jecon: %f \t, Jcontrol: %f, \t Jstate: %f \n', Jecon,Jcontrol,Jstate);
Jobj.reg    = J;
Jobj.econ   = Jecon;

end
