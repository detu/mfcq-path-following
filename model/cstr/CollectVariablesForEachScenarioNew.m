function [J,g,w0,w,lbg,ubg,lbw,ubw,t0,x0_measure] = CollectVariablesForEachScenarioNew(optProblem, N, u0, tmeasure, xmeasure, scrNo)

%COLLECTVARIABLESFOREACHSCENARIO Summary of this function goes here
% 
% [OUTPUTARGS] = COLLECTVARIABLESFOREACHSCENARIO(INPUTARGS) Explain usage here
% 
% Examples: 
% 
% Provide sample usage code here
% 
% See also: List related files here

% $Author: suwartad $	$Date: 2018/05/03 04:31:50 $	$Revision: 0.1 $
% Copyright: Process Control Group - NTNU Trondheim 2018

global pQ pC scenario;
% xf = xmeasure;
% uf = u0(1);
% mpciter = 1;  

    
%% scenario tree
% apply different price values for different realizations
switch scrNo
    case 1
        % realization 1 : -10%
        %pC    = 2;
        %pQ    = 0.5;
        %pC = pC - 0.1*pC;
        %pQ = pQ - 0.1*pQ;
        scenario = 1;
        
    case 2
        % realization 2 : -5%
        %pC    = 2;
        %pQ    = 0.5;
        %pC = pC - 0.05*pC;
        %pQ = pQ - 0.05*pQ;
        scenario = 2;
        
    case 3
        % realization 3 : 5%
        %pC    = 2;
        %pQ    = 0.5;
        %pC = pC + 0.05*pC;
        %pQ = pQ + 0.05*pQ;
        scenario = 3;
        
    case 4
        % realization 4 : 10%
        %pC    = 2;
        %pQ    = 0.5;
        %pC = pC + 0.1*pC;
        %pQ = pQ + 0.1*pQ;
        scenario = 4;
end

%% building an NLP problem for each realization
% obtain new initial value
[t0, x0] = measurementInitialValue ( tmeasure, xmeasure );

%x0         = max(min(x0,1.0),0); % restrict to boundaries
x0_measure =  x0;    % without noise
%x0_measure = max(min(x0_measure,1.0),0); % restrict to boundaries


% ideal NMPC:
%[primalNLP, ~, lb, ub, ~, params, elapsedtime] = solveOptimalControlProblem(optProblem, N, x0, u0, x0_measure);
% seharusnya !
[J,g,w0,w,lbg,ubg,lbw,ubw] = buildOptimalControlProblem(optProblem, N, x0, u0, x0_measure);


end
