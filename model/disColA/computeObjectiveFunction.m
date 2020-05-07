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
%xf(1) = 10*xf(1);          % scaled bottom concentration
u_opt = Xinit(85:89);

load Qmax.mat;

F_0 = 0.3;
% prices
pf = 1;
pV = 0.02;
pB = 2;
pD = 0;
%Qmax(1)    = 10;

% Qmax(1:42)=1;
% Qmax(43:84)=0;
% 
% %Qmax(43)=0.11;
% 
% Qmax(83)=0.11;
% %Qmax(83)= 0.1;

%Jcontrol   = (Qmax(nx+1:nx+nu,1).*(u - u_opt))' * (u - u_opt);
%Jstate     = 0.5.*(Qmax(1:nx,1).*(x - xf))' * (x - xf);
%Jstate     = 0.5.*(x - xf)' * (x - xf);
Jstate     = (Qmax(1:nx,1).*(x - xf))' * (x - xf);
%Jstate      = 1*(x(83) - xf(83))^2;
Jecon      = (pf*F_0 + pV*u(2) - pB*u(5));
%Jecon      = (pf*F_0 + pV*u(2) - pB*u(5)*x(1));
%Jecon      = pV*u(2);
%L          = Jcontrol + Jstate + Jecon;
L          = Jecon + Jstate;
%L          = Jecon;

% Jstate     = (Qmax(1:84,1).*(x - xf))' * (x - xf);
% L          = Jstate;

end
