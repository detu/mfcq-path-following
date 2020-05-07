clear
clc

%addpath ('C:\Users\dineshk\CasADi\casadi-matlabR2014b-v3.1.0-rc1')
import casadi.*


%% Set parameters, initial values, upper and lower bounds

% All the parameter values are defined inside this function
par          = GasLift_Param;
par.GOR_case = [par.GOR-par.GOR_var,par.GOR+[-1,0;0,1]*par.GOR_var ,par.GOR,par.GOR+[1,0;0,-1]*par.GOR_var,par.GOR+par.GOR_var]; % 3 scenarios
par.nS       = 5;

% import initial conditions
[dx0,z0,u0,lbx,lbz,lbu,ubx,ubz,ubu] = GasLift_Initialization_bounds(par);

%% Simulation Setup
par.N       = 24;           % no. of samples in prediction horizon
par.T       = 2*3600;       % 2 hrs predition horizon
par.tf      = par.T/par.N;  % each sample is 5 min
par.tSim    = 1;            % Shorter simulation time for EKF
par.nIter   = 1*60;         % number of closed loop iterations (60-->5h of CL simulation)
nIter = par.nIter;

PlotResult     = 0;        % plot prediction horizon at each sampling interval
EKF            = 0;        % use EKF? (1=Yes; 0=No)

[f,w_gl_0,~]   = NominalModel(par);     % steady-state optimization is here !
[F,w_gl_SP,~]  = WellSimulator(par);  % why there is another optimization here ?
yIndex         = [1:4,7,8,15:20];

disp(' ------------------------------------------- ')
disp(' ------------------------------------------- ')
disp([' Nominal Steady state optimum'])
disp(' ')
disp(['Well 1: ' num2str(w_gl_0(1)) ' [kg/s]'])
disp(['Well 2: ' num2str(w_gl_0(2)) ' [kg/s]'])
disp(['Total : ' num2str(sum(w_gl_0)) ' [kg/s]'])
disp(' ------------------------------------------- ')
disp([' True Steady state optimum'])
disp(' ')
disp(['Well 1: ' num2str(w_gl_SP(1)) ' [kg/s]'])
disp(['Well 2: ' num2str(w_gl_SP(2)) ' [kg/s]'])
disp(['Total : ' num2str(sum(w_gl_SP)) ' [kg/s]'])
disp(' ------------------------------------------- ')
disp(' ------------------------------------------- ')

par.dx0 = dx0;
par.z0  = z0;
xf      = dx0;
zf      = z0;
u_in    = u0;
%%
par.nu = par.n_w;
par.nz = 24;
par.nx = 6;

nx = par.nx;
nz = par.nz;
nu = par.nu;

[solver,par] = buildNLP5(f,par);
nS = par.nS;
w_opt = par.w0;
%%
d = 3;
P1 = 1;
P2 = 1;
P3 = 1;
P4 = 1;
P5 = 1;
ws_val = [1;1;1;1;1];
GOR_scen = par.GOR_case;

for sim_k = 1:nIter
    
    % True realization of GOR used in Simulation
    if sim_k<16
        RN = 0;
    else if sim_k>=30 && sim_k<45
            RN = 0.5;
        else if sim_k>=45
                RN = 1;
            else
                RN = 0.5 - 0.5*(30-sim_k)/15;
            end
        end
    end
    GOR_real = par.GOR + RN.*par.GOR_var;
    GOR_true(:,sim_k) = GOR_real;
    
    nlp_par = vertcat(par.GOR_case(:,1),par.GOR_case(:,2),par.GOR_case(:,3),par.GOR_case(:,4),par.GOR_case(:,5),...
                      u_in,u_in,u_in,u_in,u_in,xf,xf,xf,xf,xf,ws_val);
    
    tic;
    sol = solver('x0',w_opt,'p',nlp_par,'lbx',par.lbw,'ubx',par.ubw,'lbg',par.lbg,'ubg',par.ubg);
    t(sim_k) = toc;
    
    % Extract Solution
    w_opt = full(sol.x);
    n_w_i = nx+(nu + (nx+nz)*d+nx)*(par.N) ;
    
    for i = 1:nS
        u_opt1(:,i) = [w_opt((i-1)*n_w_i+(nx+1):nu+d*(nx+nz)+nx:i*n_w_i);NaN];
        u_opt2(:,i) = [w_opt((i-1)*n_w_i+(nx+2):nu+d*(nx+nz)+nx:i*n_w_i);NaN];
        w_po1(:,i)  = [w_opt((i-1)*n_w_i+nx+(nu+nx+19):nu+d*(nx+nz)+nx:i*n_w_i);NaN];
        w_po2(:,i)  = [w_opt((i-1)*n_w_i+nx+(nu+nx+20):nu+d*(nx+nz)+nx:i*n_w_i);NaN];
        w_pg1(:,i)  = [w_opt((i-1)*n_w_i+nx+(nu+nx+17):nu+d*(nx+nz)+nx:i*n_w_i);NaN];
        w_pg2(:,i)  = [w_opt((i-1)*n_w_i+nx+(nu+nx+18):nu+d*(nx+nz)+nx:i*n_w_i);NaN];
    end
    
    u_in_1 = u_opt1(1,1);
    u_in_2 = u_opt2(1,1);
    u_in = [u_in_1;u_in_2];
    
    flag = solver.stats();
    disp([flag.return_status ' ' num2str(sim_k)])
    
    for plant_sim = 1:par.tf/par.tSim
        Fk = F('x0',xf,'z0',zf,'p',[u_in;GOR_real]);
        
        % set new initial values for the next iteration
        xf =  full(Fk.xf) + 0.0*randn(1);
        zf = full(Fk.zf) + 0.0*randn(1);
        J_real(sim_k) = full(Fk.qf);
        
        y_meas = zf(yIndex);
        
       NMPC.ws_val(:,sim_k) = ws_val;
       NMPC.GOR_upd1(:,sim_k) = par.GOR_case(:,1);
       NMPC.GOR_upd2(:,sim_k) = par.GOR_case(:,2);
       NMPC.GOR_upd3(:,sim_k) = par.GOR_case(:,3);
       NMPC.GOR_upd4(:,sim_k) = par.GOR_case(:,4);
       NMPC.GOR_upd5(:,sim_k) = par.GOR_case(:,5);
       
    end
    u_in_MPC1(sim_k)  = u_in(1);
    u_in_MPC2(sim_k)  = u_in(2);
    NMPC.u_in_MPC = [u_in_MPC1;u_in_MPC2];
    
    w_po_MPC1(sim_k)  = full(Fk.zf(19));
    w_po_MPC2(sim_k)  = full(Fk.zf(20));
    NMPC.w_po_MPC = [w_po_MPC1;w_po_MPC2];
    
    w_pg_MPC1(sim_k)  = full(Fk.zf(17));
    w_pg_MPC2(sim_k)  = full(Fk.zf(18));
    NMPC.w_pg_MPC = [w_pg_MPC1;w_pg_MPC2];
    
    tgrid = linspace(0,par.T,par.N+1)./60;
    tpast = linspace(-(sim_k-1)*300,0,length(u_in_MPC1))./60;
    
    
    %% plot iterations
    PlotResult = 0;
    if PlotResult || rem(sim_k,1) ==0
        
        figure(56)
        clf
        subplot(211)
        hold all
        plot(tpast,w_pg_MPC1,'b')
        plot(tpast,w_pg_MPC2,'r')
        plot(tpast,w_pg_MPC1+w_pg_MPC2,'k')
        plot(tgrid,w_pg1(:,1),'b-.')
        plot(tgrid,w_pg1(:,2),'b--')
        plot(tgrid,w_pg1(:,3),'b--')
        plot(tgrid,w_pg2(:,1),'r-.')
        plot(tgrid,w_pg2(:,2),'r--')
        plot(tgrid,w_pg2(:,3),'r--')
        plot(tgrid,w_pg1(:,1)+w_pg2(:,1),'k-.')
        plot(tgrid,w_pg1(:,2)+w_pg2(:,2),'k--')
        plot(tgrid,w_pg1(:,3)+w_pg2(:,3),'k--')
        plot([tpast,tgrid],[par.QgMax.*ones(size(tpast)),par.QgMax.*ones(size(tgrid))],'m:')
        plot([0,0],[0,10],'k:')
        xlim([-1*300,120])
        grid on
        leg = legend('well 1','well 2','Total Gas rate');
        set(leg,'Interpreter','latex')
        ylabel ('Gas rate [kg/s]','Interpreter','LaTex')
        xlabel ('sample instant N','Interpreter','LaTex')
        %         title ([OptCase ' Optimization'],'Interpreter','LaTex')
        
        
        subplot(212)
        hold all
        stairs(tpast,u_in_MPC1,'b')
        stairs(tpast,u_in_MPC2,'r')
        stairs(tpast,u_in_MPC1+u_in_MPC2,'k')
        stairs(tgrid,u_opt1(:,1),'b-.')
        stairs(tgrid,u_opt1(:,2),'b--')
        stairs(tgrid,u_opt1(:,3),'b--')
        stairs(tgrid,u_opt2(:,1),'r-.')
        stairs(tgrid,u_opt2(:,2),'r--')
        stairs(tgrid,u_opt2(:,3),'r--')
        stairs([0,0],[0,7],'k:')
        xlim([-300,120])
        grid on
        leg = legend('well 1','well 2','Total GL rate');
        set(leg,'Interpreter','latex')
        ylabel ('Gaslift rate [kg/s]','Interpreter','LaTex')
        xlabel ('sample instant N','Interpreter','LaTex')
    end
end





