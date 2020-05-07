function NMPCDistCstr
%NMPCDISTCSTR Summary of this function goes here
% 
% NMPC for CSTR + Distillation Column A 
%
% [OUTPUTARGS] = NMPCDISTCSTR(INPUTARGS) Explain usage here
% 
% Examples: 
% 
% Provide sample usage code here
% 
% See also: List related files here

% $Author: suwartad $	$Date: 2016/06/30 11:42:49 $	$Revision: 0.1 $
% Copyright: Process Control Group - NTNU Trondheim 2016

import casadi.* 
format long;
%diary('run15092017-0946.txt');
%diary('runIpopt15092017.txt');
%diary('runPfNmpc-15092017-50.txt');
global N;
% number of mpc iteration
mpciterations = 150;
%mpciterations = 2;
%mpciterations = 5;
%mpciterations = 20;
%mpciterations = 30;
%mpciterations = 50;   % 50 x 3
% number of prediction horizon
%N             = 180;
%N             = 30;
%N             = 90;
%N             = 1000;
%N             = 240;
%N             = 300;
%N             = 60;
N             = 45;
%N             = 210;
%N             = 3;
%N             = 6;
%N             = 8;
%N             = 10;
%N             = 15;
%N             = 20;
%N             = 1100;
% sampling time
T             = 1;  % [minute]
%T             = 5;
%T             = 100;
% initial controls (different initial conditions)
%load Xinit40.mat;  
%load Xinit32.mat;
%load Xinit29.mat;  
%load Xinit28.mat;
%load Xinit30.mat
% load Xinit31.mat;
%load Xinit305.mat
%load Xinit301.mat
load Xinit3075.mat
%load Xinit3065.mat;
%u0            = Xinit40(85:89);
%u0            = Xinit31(85:89);
%u0            = Xinit305(85:89);
%u0            = Xinit301(85:89);
%u0            = Xinit29(85:89);
%u0            = Xinit32(85:89);   % not ok
%u0            = Xinit28(85:89);
%u0            = Xinit30(85:89);
%u0            = Xinit3075(85:89);
u0            = Xinit3075(85:89);
%u0            = Xinit3065(85:89);
u0            = repmat(u0,1,N);
%u0            = repmat(u0,1,N)./10; %scaling for controls
% get initial measurement (states) at time T = 0.
tmeasure      = 0.0;
%xmeasure      = Xinit40(1:84);
%xmeasure      = Xinit29(1:84);
%xmeasure      = Xinit30(1:84);
%xmeasure      = Xinit3075(1:84);
%xmeasure      = Xinit32(1:84);
%xmeasure      = Xinit28(1:84);
%xmeasure      = Xinit31(1:84);
%xmeasure      = Xinit305(1:84);
%xmeasure      = Xinit301(1:84);
xmeasure      = Xinit3075(1:84);
%xmeasure      = Xinit3065(1:84);

%load Xopt31.mat;
%load Xopt3175.mat;
%load XoptWithoutLevelControl.mat
%load Xopt31Ub65.mat;
% load Xopt31Ub64.mat;
% xGuess = Xopt;

load Xopt3175Horizon90.mat
xGuess = Xopt;

% load Xopt3164Horizon90.mat;
% xGuess = Xopt;

% either call iNMPC 
%[~, xmeasureAll, uAll, obj, optRes, params, runtime] = iNmpc(@optProblem, @system, mpciterations, N, T, tmeasure, xmeasure, u0);
% % % %save iNmpc.mat xmeasureAll uAll;   % without noise
% % % %save iNmpcNoise.mat xmeasureAll uAll;
% % % xmeasureAll_1pct = xmeasureAll;
% % % uAll_1pct        = uAll;
% % % save iNmpcNoise_1pct.mat xmeasureAll_1pct uAll_1pct;
% save iNmpc32.mat xmeasureAll uAll;
% %save iNmpc32-50.mat xmeasureAll uAll;
%save iNmpc32-14.mat xmeasureAll uAll obj;
%save iNmpc32-14Wo.mat xmeasureAll uAll obj;

% new iNmpcDual that collect dual variable as well 
[~, xmeasureAll, uAll, obj, primalRes, dualRes, params, runtime] = iNmpcDual(@optProblem, @system, mpciterations, N, T, tmeasure, xmeasure, u0, xGuess);
% save iNmpcDual.mat xmeasureAll uAll obj primalRes dualRes params runtime;

% compute number of active-bound constraints
%[~, xmeasureAll, uAll, obj, optRes, params, runtime, numActiveBoundRecord] = iNmpcActiveBound(@optProblem, @system, mpciterations, N, T, tmeasure, xmeasure, u0);


% or pf-NMPC
%[~, xmeasureAll_pf, uAll_pf, obj_pf, optRes_pf, params_pf, runtime_pf, etaRecord, numActiveBoundRecord, nActiveChange] = pfNmpc(@optProblem, @system, mpciterations, N, T, tmeasure, xmeasure, u0);
% % save pfNmpc.mat xmeasureAll_pf uAll_pf; % without noise 
% % save pfNmpcNoise.mat xmeasureAll_pf uAll_pf;
% % xmeasureAll_pf_1pct = xmeasureAll_pf;
% % uAll_pf_1pct = uAll_pf;
% % save pfNmpcNoise_1pct.mat xmeasureAll_pf_1pct uAll_pf_1pct;
% %save pfNmpcWoHc.mat xmeasureAll_pf uAll_pf etaRecord numActiveBoundRecord;
% %save pfNmpcWithoutHc.mat xmeasureAll_pf uAll_pf etaRecord numActiveBoundRecord;
% %save pfNmpcPC.mat xmeasureAll_pf uAll_pf etaRecord numActiveBoundRecord;
% %save pfNmpc32WithoutLP.mat xmeasureAll_pf uAll_pf etaRecord numActiveBoundRecord;
% %save pfNmpc32WithLP.mat xmeasureAll_pf uAll_pf etaRecord numActiveBoundRecord;
% save pfNmpc32.mat xmeasureAll_pf uAll_pf etaRecord numActiveBoundRecord;

%save pfNmpc32-14.mat xmeasureAll_pf uAll_pf etaRecord numActiveBoundRecord obj_pf;
%save pfNmpc32-14Wo.mat xmeasureAll_pf uAll_pf etaRecord numActiveBoundRecord obj_pf;


%[~, xmeasureAll_pf, uAll_pf, obj_pf, primalRes_pf, dualRes_pf, params_pf, runtime_pf, etaRecord, numActiveBoundRecord, nActiveChange] = pfNmpcDual(@optProblem, @system, mpciterations, N, T, tmeasure, xmeasure, u0, xGuess);
%save pfNmpcDual.mat xmeasureAll_pf uAll_pf obj_pf primalRes_pf dualRes_pf params_pf runtime_pf etaRecord numActiveBoundRecord nActiveChange;

%[~, xmeasureAll_pfm, uAll_pfm, obj_pfm, primalRes_pfm, dualRes_pfm, params_pfm, runtime_pfm, etaRecordm, numActiveBoundRecordm, nActiveChangem] = pfNmpcDual(@optProblem, @system, mpciterations, N, T, tmeasure, xmeasure, u0);
%save pfNmpcDualm.mat xmeasureAll_pfm uAll_pfm obj_pfm primalRes_pfm dualRes_pfm params_pfm runtime_pfm etaRecordm numActiveBoundRecordm nActiveChangem;

%diary('off');

%save OnlyEconomicObj.mat; 
%save OnlyEconomicObjWithNoise.mat

%save OnlyEconomicObjUb73.mat;
%save OnlyEconomicObjUb74.mat;
%save OnlyEconomicObjUb74WithNoise.mat;
%save OnlyEconomicObjUb74WithNoiseDifferentBounds.mat;
%save OnlyEconomicObj13022019.mat;
%save OnlyEconomicObj14022019WithNoise.mat
%save OnlyEconomicObj14022019WithoutNoise.mat
%save OnlyEconomicObj14022019WithoutNoiseWithoutTopLevelControlled.mat

%save OnlyEconomicWithoutNoisePF.mat
%save OnlyEconomicWithNoisePF.mat

%save OnlyEconomicObj15022019WithoutNoise.mat

%save OnlyEconomicObj17022019WithoutNoise.mat
%save OnlyEconomicObj17022019WithoutNoisePF.mat

%save OnlyEconomicObj17022019WithNoise.mat
%save OnlyEconomicObj17022019WithNoisePF.mat

%save NormalEconomicObj18022019WithoutNoise.mat
%save NormalEconomicObj18022019WithNoise.mat

%save NormalEconomicObj18022019WithNoisePF.mat

%save NormalEconomicObj19022019WithNoise.mat;
%save NormalEconomicObj19022019WithoutNoise.mat;

%save NormalEconomicObj23022019WithNoisePF.mat;
%save NormalEconomicObj24022019WithoutNoisePF.mat;
keyboard;

%% THE CODE BELOW IS JUST FOR PLOTTING
% load CstrDistXinit.mat;
% xf    = Xinit(1:84);
% u_opt = Xinit(85:89);
% 
% nu      = size(u0,1);
% uAll    = reshape(uAll,nu,mpciterations);
% uAll_pf = reshape(uAll_pf,nu,mpciterations);
% 
% % add initial control
% uAll    = [u0(:,1) uAll];
% uAll_pf = [u0(:,1) uAll_pf];
% 
% % add initial states
% xmeasureAll    = horzcat(xmeasure,xmeasureAll);
% xmeasureAll_pf = horzcat(xmeasure,xmeasureAll_pf);
% 
% %global nk;
% %x = linspace(1,nk*N*T,mpciterations);
% x = linspace(1,mpciterations,mpciterations/T);
% xi = [0 x];
% 
% close all; % close all figures
% 
% figure(1);
% clf;
% % figure('Units', 'pixels', ...
% %     'Position', [100 100 500 375]);
% hold on;
% hV_nlp = plot(x,obj,'LineWidth',6.0,'Color','g');
% hold on;
% hV_pf  = plot(x,obj_pf,'LineWidth',1.0,'Color','k');
% hTitle  = title ('Objective functions comparison PF - NLP');
% hXLabel = xlabel('Number of MPC iteration [-]'             );
% hYLabel = ylabel('Objective function [-]'                  );
% hLegend = legend( ...
%   [hV_nlp, hV_pf],  ...
%   'NLP: Obj. Func.' , ...
%   'PF:  Obj. Func.' , ...
%   'location', 'NorthWest' );
% 
% set( gca                       , ...
%     'FontName'   , 'Helvetica' );
% set([hTitle, hXLabel, hYLabel], ...
%     'FontName'   , 'AvantGarde');
% set([hLegend, gca]             , ...
%     'FontSize'   , 8           );
% set([hXLabel, hYLabel]  , ...
%     'FontSize'   , 10          );
% set( hTitle                    , ...
%     'FontSize'   , 12          , ...
%     'FontWeight' , 'bold'      );
% 
% set(gca, ...
%   'Box'         , 'off'     , ...
%   'TickDir'     , 'out'     , ...
%   'TickLength'  , [.02 .02] , ...
%   'XMinorTick'  , 'on'      , ...
%   'YMinorTick'  , 'on'      , ...
%   'YGrid'       , 'on'      , ...
%   'XColor'      , [.3 .3 .3], ...
%   'YColor'      , [.3 .3 .3], ...
%   'LineWidth'   , 1         );
% 
% set(gcf, 'PaperPositionMode', 'auto');
% % print -depsc2 ObjFunc.eps
% 
% 
% 
% figure(2);
% clf;
% plot(xi,xmeasureAll(1,:),'o','LineWidth',6.0,'Color','g');
% hold on;
% plot(xi,xmeasureAll_pf(1,:),'*','LineWidth',1.0,'Color','k');
% hold on;
% plot(xi, xf(1)*ones(1,mpciterations+1),'-r');
% title ('x(1) comparison PF - NLP');
% xlabel('Time [minute]'                      );
% ylabel('Concentration [-]'                  );
% 
% 
% figure(3);
% clf;
% %plot(xmeasureAll(2,:),'LineWidth',2.5);
% plot(xi,xmeasureAll(21,:),'LineWidth',6.0,'Color','g');
% hold on;
% plot(xi,xmeasureAll_pf(21,:),'LineWidth',1.0,'Color','k');
% hold on;
% plot(xi, xf(21)*ones(1,mpciterations+1),'-r');
% title ('x(21) comparison PF - NLP');
% xlabel('Time [minute]'                      );
% ylabel('Concentration [-]'                  );
% 
% 
% figure(4);
% clf;
% plot(xi,xmeasureAll(41,:),'LineWidth',6.0,'Color','g');
% hold on;
% plot(xi,xmeasureAll_pf(41,:),'LineWidth',1.0,'Color','k');
% hold on;
% plot(xi, xf(41)*ones(1,mpciterations+1),'-r');
% title ('x(41) comparison PF - NLP');
% xlabel('Time [minute]'                      );
% ylabel('Concentration [-]'                  );
% 
% figure(5);
% clf;
% %plot(xmeasureAll(3,:),'LineWidth',2.5);
% plot(xi,xmeasureAll(42,:),'LineWidth',6.0,'Color','g');
% hold on;
% plot(xi,xmeasureAll_pf(42,:),'LineWidth',1.0,'Color','k');
% hold on;
% plot(xi, xf(42)*ones(1,mpciterations+1),'-r');
% title ('x(42) comparison PF - NLP');
% xlabel('Time [minute]'                      );
% ylabel('Concentration [-]'                  );
% 
% figure(6);
% clf;
% %plot(xmeasureAll(3,:),'LineWidth',2.5);
% plot(xi,xmeasureAll(82,:),'LineWidth',6.0,'Color','g');
% hold on;
% plot(xi,xmeasureAll_pf(82,:),'LineWidth',1.0,'Color','k');
% hold on;
% plot(xi, xf(82)*ones(1,mpciterations+1),'-r');
% title ('x(82) comparison PF - NLP');
% xlabel('Time [minute]'                      );
% ylabel('Hold-up [-]'                        );
% 
% figure(7);
% clf;
% plot(xi,uAll(1,:),'LineWidth',6.0,'Color','g');
% hold on;
% plot(xi,uAll_pf(1,:),'LineWidth',1.0,'Color','k');
% hold on;
% plot(xi, u_opt(1,:)*ones(1,mpciterations+1),'-r');
% title ('u(1)-LT: control input comparison PF - NLP');
% xlabel('Time [minute]'                          );
% ylabel('LT [m^3/minute]'                        );
% 
% figure(8);
% clf;
% plot(xi,uAll(2,:),'LineWidth',6.0,'Color','g');
% hold on;
% plot(xi,uAll_pf(2,:),'LineWidth',1.0,'Color','k');
% hold on;
% plot(xi, u_opt(2,:)*ones(1,mpciterations+1),'-r');
% title ('u(2)-VB: control input comparison PF - NLP');
% xlabel('Time [minute]'                          );
% ylabel('VB [kmol/minute]'                        );
% 
% figure(9);
% clf;
% plot(xi,uAll(3,:),'LineWidth',6.0,'Color','g');
% hold on;
% plot(xi,uAll_pf(3,:),'LineWidth',1.0,'Color','k');
% hold on;
% plot(xi, u_opt(3,:)*ones(1,mpciterations+1),'-r');
% title ('u(3)-F: control input comparison PF - NLP');
% xlabel('Time [minute]'                          );
% ylabel('F [kmol/minute]'                        );
% 
% figure(10);
% clf;
% plot(xi,uAll(4,:),'LineWidth',6.0,'Color','g');
% hold on;
% plot(xi,uAll_pf(4,:),'LineWidth',1.0,'Color','k');
% hold on;
% plot(xi, u_opt(4,:)*ones(1,mpciterations+1),'-r');
% title ('u(4)-D: control input comparison PF - NLP');
% xlabel('Time [minute]'                          );
% ylabel('D [kmol/minute]'                        );
% 
% figure(11);
% clf;
% plot(xi,uAll(5,:),'LineWidth',6.0,'Color','g');
% hold on;
% plot(xi,uAll_pf(5,:),'LineWidth',1.0,'Color','k');
% hold on;
% plot(xi, u_opt(5,:)*ones(1,mpciterations+1),'-r');
% title ('u(5)-B: control input comparison PF - NLP');
% xlabel('Time [minute]'                          );
% ylabel('B [kmol/minute]'                        );
% 
% keyboard;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Definition of the NMPC functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = system(t, x, u, T)
    
    global uc;
    uc = u;
    [~,x_out] = ode15s('cola_lv_cstr',[t t+T], x);
    lengthx   = size(x_out); 
    y         = x_out(lengthx(1),:)'; 
    
end

function [J,g,w0,w,lbg,ubg,lbw,ubw,params] = optProblem(x, u, N, x0_measure)    %add prediction horizon 
    import casadi.*
    
    %xc = 10;
    xc = 1;
    %scaling for bottom concentration
    %x(:,1) = xc/10;
    
    % the model
    NT = 41;
    global Uf;
    Uf = 0.31;           % Feeding rate F_0
    %Uf = 0.30;
    %Uf = 0.3101; % OK
    %Uf = 0.3105; % OK
    
    % invoke the model
    [~,state,xdot,inputs] = DistColACstr(Uf);
    %[~,state,xdot,inputs] = DistColACstrScaled(Uf);  % with scaled equations
    %f = Function('f',{state,inputs}, {xdot});
    
    % objective function value
    obj = computeObjectiveFunction(state,inputs);
 
    f = Function('f',{state,inputs}, {xdot,obj});
    
    % bound constraints
    VB_max = 4.008;
    %xB_max = 0.1;
    xB_max = xc/10;    %scaled bottom concentration
    
    
    % State bounds and initial guess
    x_min =  zeros(84,1);  % try without epsilon here, later put epsilon
    %x_min =  1e-6*ones(84,1);  % try without epsilon here, later put epsilon
    x_max =  ones(84,1);
    
%     %x_min(84) = 0.3;
%     x_min(84) = 0.74;
     x_max(1)  = xB_max; % scaled bottom concentration
%     %x_max(1)  = 0.1;
%     %x_min(1)  = 0.09999;
%     %x_max(84) = 0.7;
%     %x_max(84) = 0.74;
%     %x_min(84) = 0.7;
%     x_min(84) = 0.7;  % with measurement noise
     x_max(84) = 0.75;
%     %x_max(84) = 0.72;
%     
% %     x_min(43) = 0.5;
% %     x_max(43) = 0.5;
%     %x_min(43) = 0.3;
%     %x_min(43) = 0.5;
%     x_min(43) = 0.4; % with measurement noise
%     x_max(43) = 0.6;
%     x_min(43) = 0.45; % with measurement noise
%     x_max(43) = 0.55;
%     
% %     x_min(83) = 0.5;
% %     x_max(83) = 0.6;
% 
% %     x_min(83) = 0.5;
% %     x_max(83) = 0.5;
    
    % Control bounds
    u_min = [0.1; 0.1; 0.1; 0.1; 0.1];
    u_max = [10; VB_max; 10; 1.0; 1.0];
    
    % soft constraint bounds
    sc_min = 0;
    sc_max = 1e6;
    
    % compact bound variable for ease of function invocation 
    params.bound.x_min  = x_min;
    params.bound.x_max  = x_max;
    params.bound.u_min  = u_min;
    params.bound.u_max  = u_max;
    params.bound.sc_min = sc_min;
    params.bound.sc_max = sc_max;
    
    % Construct objective function
    load CstrDistXinit.mat;
    xf    = Xinit(1:84);
    xf(1) = xc*xf(1);          % scaled bottom concentration
    %xf    = Xinit(1:84)./0.1; %scaled state
    u_opt = Xinit(85:89);
    %u_opt = Xinit(85:89)./10; %scaled control
    
    % prices
    pf = 1; 
    pV = 0.02;
    pB = 2;
    pD = 0;

    % compact price variable
    params.price.pf = pf;
    params.price.pV = pV;
    params.price.pB = pB;
    params.price.pD = pD;
    params.price.F_0= Uf;
    
    % dimensions
    global nx nu nk d tf ns;
    nx = 84;   % CSTR + Distillation Column A
    nu = 5;    % LT, VB, F, D, B
    nk = 1;
    tf = 1;   % in [minutes]
    %tf = 5;
    %tf = 100;
    h  = tf/nk;
    ns = 0;
    
    % compact model variable
    params.model.NT = NT;
    params.model.f  = f;
    params.model.xdot_val_rf_ss = xf;
    params.model.x  = x;
    params.model.u_opt = u_opt;
    % duplicate u0
    params.model.u  = repmat(u,1,nk); 


    % preparing collocation matrices
    [B,C,D,d] = collocationSetup();
    
    % compact collocation variable
    params.colloc.C = C;
    params.colloc.D = D;
    params.colloc.h = h;
    params.colloc.B = B;

    % start with an empty NLP
    w   = {};      % decision variables contain both control and state variables
    w0  = [];      % initial guess
    lbw = [];      % lower bound for decision variable
    ubw = [];      % upper bound
    J   = 0;       % objective function
    g   = {};      % nonlinear constraint
    lbg = [];      % lower bound for nonlinear constraint
    ubg = [];      % upper bound

    %delta_time = 60; % [minute] convert second to minute
    delta_time = 1;
    alpha = 1;
    beta  = 1;
    gamma = 1;
    
    % compact weight variable
    params.weight.delta_time = delta_time;
    params.weight.alpha      = alpha;
    params.weight.beta       = beta;
    params.weight.gamma      = gamma;
    
    % "Lift" initial conditions
    X0  = MX.sym('X0', nx);
    w   = {w{:}, X0};
    lbw = [lbw; x_min];
    ubw = [ubw; x_max];
    
    %w0  = [w0; x(1,1:nx)'];
    w0  = [w0; x0_measure];
    
    g   = {g{:}, X0 - x0_measure};  % USE MEASUREMENT HERE !
    lbg = [lbg; zeros(nx,1)];
    ubg = [ubg; zeros(nx,1)];

    % formulate the NLP
    Xk = X0;

    load Qmax.mat;
    params.Qmax = Qmax;
    
    %Qmax(1) = 0;
    %params.Qmax = 1e6*Qmax;  %WORKS!
    %params.Qmax = 1e2*Qmax;
    %params.Qmax = 5e1*Qmax;
    
    %params.Qmax = 20*Qmax;
    
    %params.Qmax = 5*Qmax;
    %params.Qmax = 10*Qmax;
    
%     % Fixed step Runge-Kutta 4 integrator
%     M  = 4; % RK4 steps per interval
%     DT = tf/nk/M;
%     Xo = MX.sym('Xo', nx);
%     U = MX.sym('U', nu);
%     X = Xo;
%     for j=1:M
%         k1 = f(X, U);
%         k2 = f(X + DT/2 * k1, U);
%         k3 = f(X + DT/2 * k2, U);
%         k4 = f(X + DT * k3, U);
%         X  = X+DT/6*(k1 +2*k2 +2*k3 +k4);
%         %Q = Q + DT/6*(k1_q + 2*k2_q + 2*k3_q + k4_q);
%     end
%     F = Function('F', {Xo, U}, {X}, {'x0','p'}, {'xf'});
    
%     % CVODES from the SUNDIALS suite
%     dae = struct('x',state,'p',inputs,'ode',xdot);
%     opts = struct('tf',h);
%     F = integrator('F', 'cvodes', dae, opts);
%     
%     params.F = F;
    
    % HERE SHOULD BE LOOP N-TIMES ACCORDING TO THE NUMBER OF PREDICTION HORIZON
    %count  = 2; % counter for state variable as initial guess
    count = 1;
    ssoftc = 0;
    for i=1:N
        [J,g,w0,w,lbg,ubg,lbw,ubw,Xk,params,count,ssoftc] = iterateOnPredictionHorizon(Xk, w, w0, lbw, ubw, lbg, ubg, g, J, params, i, count,ssoftc);
    end

    
end


function [J,g,w0,w,lbg,ubg,lbw,ubw,Xk,params,count,ssoftc] = iterateOnPredictionHorizon(Xk, w, w0, lbw, ubw, lbg, ubg, g, J, params, iter, count, ssoftc)

   import casadi.*
   global N;
   % extract compact variables
   x_min  = params.bound.x_min;
   x_max  = params.bound.x_max;
   u_min  = params.bound.u_min;
   u_max  = params.bound.u_max;
   sc_min = params.bound.sc_min;
   sc_max = params.bound.sc_max;
   
   
   f  = params.model.f;
   xdot_val_rf_ss = params.model.xdot_val_rf_ss;
   x = params.model.x;
   u = params.model.u;
   u_opt = params.model.u_opt;
   
   pf = params.price.pf;
   pV = params.price.pV;
   pB = params.price.pB;
   pD = params.price.pD;
   F_0= params.price.F_0;
   
   C = params.colloc.C;
   D = params.colloc.D;
   h = params.colloc.h;
   B = params.colloc.B;
   
   delta_time = params.weight.delta_time;
   Qmax = params.Qmax;
   
   global nx nu nk d tf ns;
   
%    F = params.F;
   

% with collocation method
   for k=0:nk-1
        % New NLP variable for the control
        Uk     = MX.sym(['U_' num2str((iter-1)*nk+k)], nu);
        w      = {w{:}, Uk};
        lbw    = [lbw; u_min];
        ubw    = [ubw; u_max];
        indexU = (iter-1)*nk + (k+1);
        w0     = [w0;  u(:,indexU)];
        
        Jcontrol   = (Qmax(nx+1:nx+nu,1).*(Uk - u_opt))' * (Uk - u_opt);
      
        % State at collocation points
        Xkj   = {};
        Skj   = {};
        Ekj   = {};
        Jcoll = 0;
        sumE  = {};
        for j=1:d
            Xkj{j} = MX.sym(['X_' num2str((iter-1)*nk+k) '_' num2str(j)], nx);
%             % soft-constraint variables
%             Skj{j} = MX.sym(['S_' num2str((iter-1)*nk+k) '_' num2str(j)], 1);
%             Ekj{j} = MX.sym(['E_' num2str((iter-1)*nk+k) '_' num2str(j)], 1);
%             % append soft-constraint variables
            w      = {w{:}, Xkj{j}};
%             w      = {w{:}, Skj{j}};
%             w      = {w{:}, Ekj{j}};

%             if j == 1
%                 % equality constraints only at the beginning finite element
%                 % for bottom level, top level, and reactor holdup
%                 xub = x_max;
%                 xlb = x_min;
%                 % bottom level
%                 xub(43) = 0.5;
%                 xlb(43) = 0.5;
%                 % top level
%                 xub(83) = 0.5;
%                 xlb(83) = 0.5;
% %                 % reactor holdup
% %                 xub(84) = 0.75;
% %                 xlb(84) = 0.74;
%                 lbw    = [lbw; xlb];
%                 ubw    = [ubw; xub];
%             else
%                 lbw    = [lbw; x_min];
%                 ubw    = [ubw; x_max];
%             end
            lbw    = [lbw; x_min];
            ubw    = [ubw; x_max];
            
%             lbw    = [lbw; x_min; sc_min; sc_min];
%             ubw    = [ubw; x_max; sc_max; sc_max];
            
            %w0     = [w0; x(iter+1,:)'];
            w0     = [w0; x(:,count)];
            
            
            %w0     = [w0; x(iter+1,:)'; 0; 0];
            count  = count + 1;
            
%             % Append soft constraint
%             g   = {g{:}, Xkj{j}(1) + Skj{j} - Ekj{j} - 1};
%             lbg = [lbg; 0];
%             %ubg = [ubg; x_max(1)];
%             ubg = [ubg; 0];
            
            Jcoll = Jcoll + (Qmax(1:nx,1).*(Xkj{j} - xdot_val_rf_ss))' * (Xkj{j} - xdot_val_rf_ss) * delta_time;
            
            %sumE = {sumE{:} Ekj{j}};
            
        end

        % Loop over collocation points
        Xk_end = D(1)*Xk;

        for j=1:d
           % Expression for the state derivative at the collocation point
           xp = C(1,j+1)*Xk;
           for r=1:d
               xp = xp + C(r+1,j+1)*Xkj{r};
           end

           % Append collocation equations
           %fj  = f(Xkj{j},Uk);
           [fj, qj] = f(Xkj{j},Uk);
           g   = {g{:}, h*fj - xp};
           lbg = [lbg; zeros(nx,1)];
           ubg = [ubg; zeros(nx,1)];

           % Add contribution to the end state
           Xk_end = Xk_end + D(j+1)*Xkj{j};
           
           J = J + B(j+1)*qj*h;
        end    

        % New NLP variable for state at end of interval
        Xk  = MX.sym(['X_' num2str((iter-1)*nk+k)], nx);
%         % soft-constraint variables
%         Sk  = MX.sym(['S_' num2str((iter-1)*nk+k) '_' num2str(j)], 1);
%         Ek  = MX.sym(['E_' num2str((iter-1)*nk+k) '_' num2str(j)], 1);
        
        w   = {w{:}, Xk};
%         w   = {w{:}, Xk, Sk, Ek};
        
        % add soft-constraint for bottom composition
        lbw = [lbw; x_min];
        ubw = [ubw; x_max];
%         lbw    = [lbw; x_min; sc_min; sc_min];
%         ubw    = [ubw; x_max; sc_max; sc_max];
        
        %w0  = [w0; x(iter+1,:)'];
        w0   = [w0; x(:,count)];
        
        %w0     = [w0; x(iter+1,:)'; 0; 0];
        count  = count + 1;

        % Add equality constraint
        g   = {g{:}, (Xk_end-Xk)};
        lbg = [lbg; zeros(nx,1)];
        ubg = [ubg; zeros(nx,1)];
        
%         % Append soft constraint
%         g   = {g{:}, Xk(1) + Sk - Ek - 1};
%         lbg = [lbg; 0];
%         %ubg = [ubg; x_max(1)];
%         ubg = [ubg; 0];
               
        Jecon  = (pf*F_0 + pV*Uk(2) - pB*Uk(5) - pD*Uk(4)) * delta_time;
        Jstate =(Qmax(1:nx,1).*(Xk - xdot_val_rf_ss))' * (Xk - xdot_val_rf_ss) * delta_time;
        
        alpha  = 1;
        beta   = 1;
        gamma  = 1;

        %J = J + alpha*Jcontrol + gamma*Jstate + beta*Jecon + Jcoll;
        
%         % add soft-constraint in the objective function
%         sumE = {sumE{:} Ek};
%         %J = J + 1e3*(sumE*ones(4,1));
%         for i=1:4
%             J = J + 1e3*sumE{i};
%         end
    end

% % with RK integrator
% for k=0:nk-1
%     %New NLP variable for the control
%     Uk     = MX.sym(['U_' num2str((iter-1)*nk+k)], nu);
%     w      = {w{:}, Uk};
%     lbw    = [lbw; u_min];
%     ubw    = [ubw; u_max];
%     indexU = (iter-1)*nk + (k+1);
%     w0     = [w0;  u(:,indexU)];
%     
%     Jcontrol   = (Qmax(nx+1:nx+nu,1).*(Uk - u_opt))' * (Uk - u_opt);
%     
%     % Integrate till the end of the interval
%     Fk = F('x0', Xk, 'p', Uk);
%     Xk_end = Fk.xf;
%     %J=J+Fk.qf;
% 
%     % New NLP variable for state at end of interval
%     Xk  = MX.sym(['X_' num2str(k+1)], nx);
%     w   = [w, {Xk}];
%     lbw = [lbw; x_min];
%     ubw = [ubw; x_max];
%     w0     = [w0; x(iter+1,:)'];
%     count  = count + 1;
% 
%     % Add equality constraint
% %     g = [g, {Xk_end-Xk}];
% %     lbg = [lbg; 0; 0];
% %     ubg = [ubg; 0; 0];
%     g   = {g{:}, (Xk_end-Xk)};
%     lbg = [lbg; zeros(nx,1)];
%     ubg = [ubg; zeros(nx,1)];
%     
%     Jecon  = (pf*F_0 + pV*Uk(2) - pB*Uk(5) - pD*Uk(4)) * delta_time;
%     Jstate =(Qmax(1:nx,1).*(Xk - xdot_val_rf_ss))' * (Xk - xdot_val_rf_ss) * delta_time;
%     
%     alpha  = 1;
%     beta   = 1;
%     gamma  = 1;
%     
%     J = J + alpha*Jcontrol + gamma*Jstate + beta*Jecon;
% end
end


