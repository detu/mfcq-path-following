% a script to plot lam_x (Lagrange multiplier for bound constraint)

load ycell.mat;
load bounds.mat;
numCell = size(ycell,2);

% extract data from ycell
for i=1:numCell
    lamX(:,i) = ycell{i}.lam_x;
end

% reshape data according to collocation point
[dual_u, dual_x] = plotStatesN(lamX(:,4), lb, ub, N);

% plot data