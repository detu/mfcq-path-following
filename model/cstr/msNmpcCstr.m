function [u_nlp_opt,t0,x0_measure,elapsednlp] = msNmpcCstr(optProblem, N, u0, tmeasure, xmeasure, scr)
%MSNMPCCSTR Summary of this function goes here
% 
% [OUTPUTARGS] = MSNMPCCSTR(INPUTARGS) Explain usage here
% 
% Examples: 
% 
% Provide sample usage code here
% 
% See also: List related files here

% $Author: suwartad $	$Date: 2018/05/14 15:00:32 $	$Revision: 0.1 $
% Copyright: Process Control Group - NTNU Trondheim 2018

import casadi.*

%% ITERATE ON THE SCENARIO TREES
% prepare optimization variables
J   = 0;
g   = {};
w0  = [];
w   = {};
lbg = [];
ubg = [];
lbw = [];
ubw = [];
nac = [];
for i=1:scr.numScr
    %for i=1:4
    % call NMPC for each scenario
    % divide into 2 realizations: +- 10%
    
    %% GATHER ALL VARIABLES FROM EACH SCENARIO
    % sNmpcCstr(optProblem, system, mpciterations, N, T, tmeasure, xmeasure, u0, scr, varargin)
    %[Tall, xmeasureAll, uAll, ObjVal, primalNLP, params, runtime] = CollectVariablesForEachScenario(optProblem, system, mpciterations, N, T, u0, tmeasure, xmeasure, i);
    [Ji,gi,w0i,wi,lbgi,ubgi,lbwi,ubwi,t0,x0_measure] = CollectVariablesForEachScenario(optProblem, N, u0, tmeasure, xmeasure, i);
    
    %% collection optimization variables for all realization
    J   = J + Ji;
    %g   = [g;gi];
    g   = {g{:},gi{:}};
    w0  = [w0;w0i];
    w   = {w{:},wi{:}};
    lbg = [lbg;lbgi];
    ubg = [ubg;ubgi];
    lbw = [lbw;lbwi];
    ubw = [ubw;ubwi];
    nac = [nac;wi{3}];  %control input
    %nac = [nac;wi{2}];   % remember control input index!
end

%% BUILD NON-ANTICIPATIVITY CONSTRAINTS (NAC)
% ingat ada steady-state constraint JANGAN DIPAKE sebagai NAC!
% let try robust horizon = 1
numNac = size(nac,1);
for i=1:(numNac-1)
    g   = {g{:}, nac(i+1,1)-nac(i,1)};
    lbg = [lbg; 0];
    ubg = [ubg; 0];
end


%% SOLVE NLP FOR ALL REALIZATIONS
prob = struct('f', J, 'x', vertcat(w{:}), 'g', vertcat(g{:}));
options = struct;
options.ipopt.tol                = 1e-12;
options.ipopt.constr_viol_tol    = 1e-12;
solver = nlpsol('solver', 'ipopt', prob, options);

% Solve the NLP
startnlp   = tic;
sol        = solver('x0', w0, 'lbx', lbw, 'ubx', ubw, 'lbg', lbg, 'ubg', ubg);
elapsednlp = toc(startnlp);
fprintf('IPOPT solver runtime = %f\n',elapsednlp);
success = strcmp(solver.stats.return_status,'Infeasible_Problem_Detected');
if (success)
    keyboard;
end

u      = full(sol.x);
lamda  = full(sol.lam_g);
objVal = full(sol.f);

%% INJECT TO PLANT
% reshape optimized variables
numCol = size(u,1)/scr.numScr;
u      = reshape(u,numCol,scr.numScr);
% loop for all scenarios for plotting purposes
u_nlp_opt  = zeros(scr.numScr,50); % 150 is hardcore
x_nlp_opt_a = zeros(scr.numScr,201); % idem
x_nlp_opt_b = zeros(scr.numScr,201);
for i= 1:scr.numScr
    uopt   = u(:,i);
    [u_nlp_opt(i,:), x_nlp_opt, ~] = arrangeOptResults(uopt, lbw(1:numCol,1), ubw(1:numCol,1), N);
    x_nlp_opt_a(i,:) = x_nlp_opt(1,:);
    x_nlp_opt_b(i,:) = x_nlp_opt(2,:);
end

end
