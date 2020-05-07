function xprime=cstr_odeWTC(t,X) 
% ODE driver for CSTR 

% global u NT %Make perturbed inputs/disturbances available to model
global uc;

% Store all inputs and disturbances
u_all = uc;

%xprime=colamod_cstr(t,X,u_all);
xprime=mod_cstrWTC(t,X,u_all);