function par = GasLift_Param()
par.n_w = 2;


par.L_w = [1500;1500];
par.H_w = [1000;1000];
par.D_w = [0.121;0.121];

par.L_bh = [500;500];
par.H_bh = [100;100];
par.D_bh = [0.121;0.121];

par.L_a = par.L_w;
par.H_a = par.H_w;
par.D_a = [0.189;0.189];

par.rho_o = [9;8].*1e2;
par.C_iv = [0.1e-3;0.1e-3];
par.C_pc = [1e-3;1e-3];

par.GOR = [0.1;0.12];
par.p_m = [20;20];
par.p_res = [150;155];
par.PI = [2.2;2.2];
par.T_a = [28;28];
par.T_w = [32;32];

par.GOR_var = [0.05;0.01];
par.rho_var = [150;25].*0;

par.R = 8.314;
par.Mw = 20e-3;
R = par.R;
Mw = par.Mw;


par.qGLMax = 50;
par.QgMax = 8;
par.qGLROC = [0.2;0.2];

qGLMax = par.qGLMax;
QgMax = par.QgMax;
qGLROC = par.qGLROC;

NyVar = 0; % measurement noise variance  
NxVar = 0; 
NzVar = 0; 

nIter =15;


useIdas = 1;
useSimulink = abs(1-useIdas);


