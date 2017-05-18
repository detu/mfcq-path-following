function [x_init, y_init, elapsedqp] = pf_pc_mfcq(problem, p_init, p_final, x_init, y_init, delta_t, lb_init, ub_init, verbose_level, N)
%PF_PC_MFCQ Summary of this function goes here
% 
% [OUTPUTARGS] = PF_PC_MFCQ(INPUTARGS) Explain usage here
%
% Implementatiion of Predictor-Corrector QP for MFCQ case
% 
% Examples: 
% 
% Provide sample usage code here
% 
% See also: List related files here

% $Author: suwartad $	$Date: 2017/02/25 21:25:12 $	$Revision: 0.1 $
% Copyright: Process Control Group - NTNU Trondheim 2017


% TO DO:
% 1. implement predictor-corrector QP with eta as stopping criteria
% 2. add LP solver 

clear prob;
sym pp;
p       = p_init;
theprob = @(pp)problem(pp);
prob    = theprob(p);
t       = 0;
alpha_2 = 0.5;
% number of iteration
iter = 0;
elapsedqp = 0;
numX    = size(x_init,1);
x0      = zeros(numX,1);   
if (verbose_level)
    fprintf('Solving problem %s \n',prob.name);
    fprintf('iteration  delta_t        t        Success\n');
end

p_0 = p_init;

% compute initial Eta
global flagDt;
flagDt = 0;
[~,g,~,~,cin,~,~,Jeq,~,~,~] = prob.obj(x_init,y_init,p, N);    % obtain derivatives information
[oldEta, ~]                 = computeEta(Jeq, g, y_init, cin);

while (t <= 1)
    
    % calculate step s
    tk   = t + delta_t;
    p_t  = (1 - tk)*p_0 + tk*p_final;
    step = p_t - p_init;
    
    % update bound constraint
    if(~isempty(lb_init))
        lb = lb_init - x_init; % (0 - 0.0999)
        ub = ub_init - x_init; %(0.1 - 0.0999)
    else
        lb = [];
        ub = [];
    end

    % solve MFCQ predictor-corrector 
    [x_init, y_init, qp_run, deltaT, newEta] = solveThreeSteps(prob, p_init, x_init, y_init, step, lb, ub, N, x0, delta_t, oldEta);  % supply initial guess
    elapsedqp = elapsedqp + qp_run;
    
    % update t
    t = t + deltaT;
    
    % update p_init (success or fail ?)
    
    % update Eta
    
    
%     if (qp_exit < 0) % QP is infeasible
%         
%         % shorten step
%         delta_t = alpha_2 * t;
%         
%         % print out iteration and FAIL
%         iter    = iter + 1;
%         success = 0;
%         if (verbose_level)
%             fprintf('%6g     %3.1e     %3.1e   %4g \n',iter, delta_t , t, success);
%         end
%         
%     else
%         % QP is feasible
%         
%         % update states, multipliers, parameter, and time step
%         x_init       = x_init + y; 
%         y_init.lam_x = y_init.lam_x - lamda.lam_x; % - (minus) because TOMLAB has different sign compared to IPOPT
% 
%         t       = t + delta_t;
%         p_init  = p_t;
%         
%         
%         % FIX STEPLENGTH
%         fprintf('delta_t: %f\n', delta_t);
%         fprintf('t: %f\n', t);
%         
%         
%         % print out iteration and SUCCESS
%         fprintf('--------------------------------------------------- \n');
%         iter    = iter + 1;
%         fprintf('iteration number: %d\n', iter);
%         success = 1;
%         if (verbose_level)
%             fprintf('%6g     %3.1e     %3.1e   %4g \n',iter, delta_t , t, success);
%         end
%         
%     end
%     if ( (1-t) <= 1e-5 )
%         break;
%     end
    
end

end
