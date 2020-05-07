function [dx0,z0,u0,lbx,lbz,lbu,ubx,ubz,ubu] = GasLift_Initialization_bounds(par)

n_w = par.n_w;
load('GasLiftSystem_forwardSimulationData.mat');

m_ga0   = [(MT1.data(end,1)).*1e-3;(MT2.data(end,1)).*1e-3];
m_gt0   = [(MT1.data(end,2)).*1e-3;(MT2.data(end,2)).*1e-3];
m_ot0   = [(MT1.data(end,3)).*1e-3;(MT2.data(end,3)).*1e-3];

p_ai0   = [(PT1.data(end,3));(PT2.data(end,3))];
p_wh0   = [(PT1.data(end,4));(PT2.data(end,4))];
p_wi0   = [(PT1.data(end,2));(PT2.data(end,2))];
p_bh0   = [(PT1.data(end,1));(PT2.data(end,1))];
rho_ai0 = [(DT1.data(end,2)).*1e-2;(DT2.data(end,2)).*1e-2];
rho_m0  = [(DT1.data(end,1)).*1e-2;(DT2.data(end,1)).*1e-2];
w_iv0   = [(FT1.data(end,4));(FT2.data(end,4))];
w_pc0   = [(FT1.data(end,3));(FT2.data(end,3))];
w_pg0   = [(FT1.data(end,1));(FT2.data(end,1))];
w_po0   = [(FT1.data(end,2));(FT2.data(end,2))];
w_ro0   = [(FT1.data(end,5));(FT2.data(end,5))];
w_rg0   = [(FT1.data(end,6)).*1e1;(FT2.data(end,6)).*1e1];

w_gl0 = [1;1];


dx0 = vertcat(m_ga0,m_gt0,m_ot0);
z0 = vertcat(p_ai0,p_wh0,p_wi0,p_bh0,rho_ai0,rho_m0,w_iv0,w_pc0,w_pg0,w_po0,w_ro0,w_rg0);
u0 = w_gl0;


m_ga_lb = [10.*1e-3;10.*1e-3];
m_gt_lb = [10.*1e-3;10.*1e-3];
m_ot_lb = [10.*1e-3;10.*1e-3];

p_ai_lb = 0.1.*ones(n_w,1);
p_wh_lb = 0.1.*ones(n_w,1);
p_wi_lb = 0.1.*ones(n_w,1);
p_bh_lb = 30.*ones(n_w,1);
rho_ai_lb = 0.01.*ones(n_w,1);
rho_m_lb = 0.01.*ones(n_w,1);
w_iv_lb = 1e-2.*ones(n_w,1);
w_pc_lb = 1e-2.*ones(n_w,1);
w_pg_lb = 1e-2.*ones(n_w,1);
w_po_lb = 1e-2.*ones(n_w,1);
w_ro_lb = 1e-2.*ones(n_w,1);
w_rg_lb = 1e-2.*ones(n_w,1);

w_gl_lb = 1e-2.*ones(n_w,1);

lbx = vertcat(m_ga_lb,m_gt_lb,m_ot_lb);
lbz = vertcat(p_ai_lb,p_wh_lb,p_wi_lb,p_bh_lb,rho_ai_lb,rho_m_lb,w_iv_lb,w_pc_lb,w_pg_lb,w_po_lb,w_ro_lb,w_rg_lb);
lbu = w_gl_lb;

m_ga_ub = 10e7.*ones(n_w,1);
m_gt_ub = 10e7.*ones(n_w,1);
m_ot_ub = 10e7.*ones(n_w,1);

p_ai_ub = 150e4.*ones(n_w,1);
p_wh_ub = 70e4.*ones(n_w,1);
p_wi_ub = 150e4.*ones(n_w,1);
p_bh_ub = 150e4.*ones(n_w,1);
rho_ai_ub = 900e4.*ones(n_w,1);
rho_m_ub = 900e4.*ones(n_w,1);
w_iv_ub = 50e4.*ones(n_w,1);
w_pc_ub = 50e4.*ones(n_w,1);
w_pg_ub = 50e4.*ones(n_w,1);
w_po_ub = 50e4.*ones(n_w,1);
w_ro_ub = 50e4.*ones(n_w,1);
w_rg_ub = 50e4.*ones(n_w,1);
w_gl_ub = 50e4.*ones(n_w,1);

ubx = vertcat(m_ga_ub,m_gt_ub,m_ot_ub);
ubz = vertcat(p_ai_ub,p_wh_ub,p_wi_ub,p_bh_ub,rho_ai_ub,rho_m_ub,w_iv_ub,w_pc_ub,w_pg_ub,w_po_ub,w_ro_ub,w_rg_ub);
ubu = w_gl_ub;
