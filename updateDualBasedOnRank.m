function [y, flag, prevRankDef, count] = updateDualBasedOnRank(y, Jeq, prevRankDef, count, oldEta)
% active bound constraints
boundMultiplier    = y.lam_x;
positiveBoundMult  = abs(boundMultiplier);
%activeIndex        = find(positiveBoundMult>1e-1);  % set active bound constraint.
%activeIndex        = find(positiveBoundMult>1e-4);
activeBoundTol     = oldEta^(0.7);
activeIndex        = find(positiveBoundMult>activeBoundTol);  % set active bound constraint.
%paramIndex         = find(activeIndex <= 84);       % remove the first 84 constraints (for distillation case!)
%activeIndex(paramIndex) = [];
numActiveBoundCons = numel(activeIndex);

[numJ,m] = size(Jeq);
JAbc     = zeros(numActiveBoundCons,m);
for i=1:numActiveBoundCons
    row         = activeIndex(i);
    JAbc(i,row) = 1;
end

Jm      = [Jeq;JAbc];


flag = 1;
nAct = numJ + numActiveBoundCons;
rk   = sprank(Jm);
%rk   = rank(full(Jm));
if rk == nAct
    flag = 0;
    return
end

rankDeficient = nAct - rk;
fprintf('rank deficient = %d\n', rankDeficient);
%if rankDeficient <= 3 % OK
%if (rankDeficient <= 2)
if (rankDeficient < 1)
    flag        = 0;
    return
end
if (prevRankDef == rankDeficient)
    count = count + 1;
else
    count = 0;
end
if (count == 4)
    flag = 0;
    return 
end
prevRankDef = rankDeficient;
% NullJ = null(full(Jm'));
% %NullJ = spnull(Jm');
% %NullJ = nullSparse(Jm');
% if isempty(NullJ)
%     keyboard;
% end
% wnJ   = NullJ(:,1);

% change with an LP problem to find a search direction 
[numX,numY] = size(Jeq);
[numA,numB] = size(JAbc);
% bound constraints
lb       = [-Inf*ones(numX+numA,1); zeros(numX+numA,1)];
ub       = [Inf*ones(numX+numA,1) ; ones(numX+numA,1)];

% objective function
f        = [zeros(numX+numA,1);-ones(numX+numA,1)];

% equality constraint
Aeq = [-[Jeq;JAbc]' [eye(numX,numY);eye(numA,numB)]'];
beq = zeros(numY,1);

%  inequality constraint
% A  = [[zeros(numX,numY);JAbc]' [zeros(numX,numY);JAbc]'];
% b  = zeros(numY,1);
A = [];
b = [];

option          = cplexoptimset;
%option.Display  = 'iter';
option.Display  = 'none';
option.lpmethod = 1; %primal simplex

[lpSol,fval,exitflag] = cplexlp(f, A, b, Aeq, beq, lb, ub, [], option );

wnJ                  = lpSol(1:numX+numA);
wnJ(find(wnJ>0))     = -eps;
%alphahat             = min(min(-(y.lam_x(activeIndex)./wnJ((numJ+1):end))),1);
alphahat             = 1;
yDummy               = [y.lam_g; y.lam_x(activeIndex)] + alphahat*wnJ;
y.lam_g              = yDummy(1:numJ);
temp                 = yDummy(numJ+1:end);
temp(find(temp<0))   = 0;
y.lam_x(activeIndex) = temp;

% % select some largest dual variables to meet the rank requirement
% rankDeficient        = nAct - rk;
% dualActiveBound      = y.lam_x(activeIndex);
% [temp,originalpos]   = sort(dualActiveBound, 'ascend');
% toBeRemoved          = originalpos(1:rankDeficient);
% dualActiveBound(toBeRemoved) = eps;
% y.lam_x(activeIndex)         = dualActiveBound;

end