function [x,y,result,itn,finaleta] = Pathfollowv2(problem,p1,p2,x0,y0,pickmul0,rate0)


% Generates a point for approximating the solution to the parametric NLP
%    min             f(x,p)
% subject to     cL <= c(x,p) <= cU
%
% prob(p) defines the problem. the dependence of f and c on p is assumed to
% be linear
% x0 is the solution to prob(p1), or an estimate thereof.
% If none is given, prob(p1) is solved as a
% standalone optimization problem.

if nargin==7
    rate = rate0;
else
    rate = 1.4;
end
if nargin==6
    pickmul = pickmul0;
else
    pickmul = 1;
end


options = optimset('Display','off','Algorithm','active-set');
%options = optimset('Display','off','Algorithm','interior-point-convex');
%optionslin = optimset('Display','off','Algorithm', 'active-set');
%optionslin = optimset('Display','off','Algorithm', 'dual-simplex');
optionslin = optimoptions('linprog','Algorithm','dual-simplex','Display','off');

result=0;
p=p1;
param = parampath;
clear prob;
sym pp;
theprob = @(pp)problem(pp);
prob = theprob(p);

m   = prob.m;
me  = prob.me;
eqs = 1:me;
ineqs = me+1:m;
n   = prob.n ;
cL  = prob.cl;      cU  = prob.cu;

if nargin == 3
    % Solve the base problem
    % Do this later
end

[f,g] = prob.obj(x0,p);
[c,J] = prob.cons( x0,p );

[etaold,~] = calceta(x0,y0,prob,f,g,c,J);

if etaold>.01
    %get a more accurate multiplier
    Act = setdiff(find(c<.01),eqs);
    
    [yact,~,exitflag,~,~] = quadprog([J(Act,:);J(eqs,:)]*[J(Act,:);J(eqs,:)]',(-2*g'*[J(Act,:);J(eqs,:)]')',[],[],[],[],[zeros(length(Act),1);-Inf*ones(length(eqs),1)],[Inf*ones(length(Act)+length(eqs),1)], y0(union(Act,eqs)),options);
    y0(Act) = yact(1:length(Act));
    y0(eqs) = yact(length(Act)+1:end);
    [etaold,~] = calceta(x0,y0,prob,f,g,c,J);
end



epsA = find(abs(y0) > etaold^param.gamma);
epsA = setdiff(epsA,eqs);
epsF = setdiff(ineqs,epsA);
Act = union(epsA,union(eqs,find(c< etaold^(param.gamma))));



converged = 0;
itn = 0;
%%4g %8.1e %s
fprintf('Solving problem %s with p=%3.1e using estimate at p=%3.1e\n',prob.name,p2(1),p1(1));
eqs = eqs';
t = 0;
deltat = param.deltat0;



fprintf('iteration  p        t        eta(x,y,p) deltat      sizeepsA    Success\n');
fprintf('%6g     %3.1e  %3.1e  %5.2e   %3.1e  %4g \n',itn,p(1),t,etaold,deltat,length(epsA));

x=x0;
y=y0;
%H = prob.hess(x,y);

etanew = etaold;
success = 0;
fail = 0;
etamax = etaold;
alphahat = 1;
ts = [0];
ys =[x(1)];
testrepeat = 0;
xs = [x'];
ys = [y'];
ts = [0];
Deltats = [deltat];
Etats   = [etanew];

%make y satisfy the necessary conditions
dcdp = prob.dcdp;

while ~converged
    itn = itn+1;
    success = 0;
    pred = 0;
    %ts = [ts;t+deltat];
    p0 = p2*(t)+p1*(1-t);
    clear p;
    p = p2*(t+deltat)+p1*(1-t-deltat);
    clear prob;
    prob = theprob(p0);
    m   = prob.m;
    n   = prob.n ;
    cL  = prob.cl;      cU  = prob.cu;
    [f,g] = prob.obj(x,p0);
    [c,J] = prob.cons(x,p0);
    keepgoingrank = 1;
    yold=y;
    if itn == 1243 || itn == 1265 || itn == 290
        1
    end

    while keepgoingrank
        nAct = length(eqs)+length(epsA);
        rk = rank(J([eqs;epsA],:));
        if rk == nAct
            keepgoingrank = 0;
            break;
        end
        NullJ = null(J([eqs;epsA],:)');
        wnJ = NullJ(:,1);
        wnJ(find(wnJ>0)) = -eps;
        alphahat = min(min(-(y([epsA])./wnJ((length(eqs)+1):end))),1);
        y([eqs;epsA]) = y([eqs;epsA])+alphahat*wnJ;
        epsA = setdiff(find(y>max(param.etamin,etaold^param.gamma)),eqs);
    end
    %norm(y-yold)
    
    
    
    HC = prob.chess(x);
    H = prob.hess(x,y,p);
    [etaold,~] = calceta(x,y,prob,f,g,c,J);
    xhat = x;
    yhat = y;
    

    
    
    
    
    
    if (etaold >=1)
        %fprintf = 'Initial eta too large, decreasing delta t\n');
        deltat = deltat*param.gamma1;
        success = 0;
    else
        %options = optimset('Display','off','TolFun',eps*etaold,'TolX',eps*etaold);
        delc = -[H -J([eqs;epsA],:)'; J([eqs;epsA],:) zeros(length(epsA)+length(eqs))]\[g-J'*y;c([eqs;epsA])];
        delxc = delc(1:n);
        delyc = delc(n+1:end);
        delycf = zeros(m,1);
        delycf([eqs;epsA]) = delyc;
        Jnew = J;
        for i=1:m
            Jnew(i,:) = Jnew(i,:)+(squeeze(HC(:,:,i))*delxc)';
        end
        [deltax,~,exitflag,~,lambda] = quadprog(H,zeros(n,1),-Jnew(epsF,:),c(epsF)+Jnew(epsF,:)*delxc+dcdp(epsF,:)*(deltat*(p2-p1)),[Jnew(eqs,:);Jnew(epsA,:)],[-dcdp(eqs,:)*(deltat*(p2-p1));-dcdp(epsA,:)*(deltat*(p2-p1))],[],[], x,options);
        %if exitflag < 0
        %    keyboard;
        %end
        dynew = zeros(m,1);
        dynew(epsF) = (lambda.ineqlin(1:(length(epsF))));
        dynew([eqs;epsA]) = -lambda.eqlin;
        prob = theprob(p);
        [f,g] = prob.obj(xhat+delxc+deltax,p);
        [c,J] = prob.cons(xhat +delxc+deltax,p);
        ynew = y+dynew+delycf;
        [etanew,z] = calceta(xhat+delxc+deltax,ynew,prob,f,g,c,J);
        
        
        if (etanew<=param.isopt) ||  (etanew <= etaold^rate && etanew< max(0.01,etamax)) || (etanew < etamax) %||  (fail>=param.maxfail && etanew< (etaold+etamax)/2)%^(1+param.gamma)
            success = 1;
            xold = x;
            yold = y;
            %ynew
            y = ynew;
            t = t+deltat;
            x = x+deltax+delxc;
            
            if (etanew <= etaold^1.6 || etanew<=param.isopt) && ~testrepeat
                deltat = min(1-t,deltat/param.gamma1);
            else
                deltat = min(1-t,deltat);
            end
            testrepeat = 0;
            etaold = etanew;
            epsA = find(abs(y) > max(param.etamin,etaold^param.gamma));
            epsA = setdiff(epsA,eqs);
            epsF = setdiff(ineqs,epsA);
            Act = [eqs;epsA;setdiff(setdiff(find(c< max(param.etamin,etaold^(param.gamma))),eqs),epsA)];
%            union(epsA,union(eqs,find(c< etaold^(param.gamma))));
            %epsA
            fail = 0;
            H = prob.hess(x,y,p);
            
            %get new multiplier
            %x = linprog(f,A,b,Aeq,beq,lb,ub,x0,options)

            lb = zeros(length(Act),1);
            lb(1:me) = -Inf;
            ub = Inf*ones(length(Act),1);

            if pickmul
            [ydel,~,exitflag] = linprog(prob.dcdp(Act,:)*deltat*(p2-p1), [J(Act,:)';-J(Act,:)'],[g+abs(z);abs(z)-g],[],[],lb,ub,y(Act),optionslin);
%             if exitflag~=1
%                 x = xold;
%                 y = yold;
%                 deltat = deltat*param.gamma1;
%                 continue;
%             end
            linobj = prob.dcdp(Act,:)*deltat*(p2-p1);
            if exitflag==1 && (linobj'*ydel < linobj'*y(Act)-sqrt(eps))      
              y(Act) = ydel;
              epsA = find(abs(y) > max(param.etamin,etaold^param.gamma));
              epsA = setdiff(epsA,eqs);
              epsF = setdiff(ineqs,epsA);
              Act = union(epsA,union(eqs,find(c< max(param.etamin^0.5,etaold^(param.gamma)))));             

            
              %y
              %epsA
            end
            end
            xs = [xs;x'];
            ts = [ts;t];
            ys = [ys;y'];
            Deltats = [Deltats;deltat];
            Etats   = [Etats;etanew];
            if etanew > etamax
                etamax = etanew;
            end
        else
            deltat = deltat*param.gamma1;
            if fail == 0
                testrepeat = 1;
            end
            fail = fail+1;
        end
    end
    
    
    
    
    if deltat <1/(100000*param.itnmax) && t<1
        result = 0;
        return
    end
    
    
    
    sizeact = length(epsA);
    % Print iteration, exit if necessary
    if mod(itn,10) == 0
        fprintf('iteration  p        t        eta(x,y,p) deltat      sizeepsA    Success\n');
    end
    fprintf('%6g     %3.1e  %3.1e  %5.2e   %3.1e  %4g        %4g\n',itn,p(1),t,etaold,deltat,sizeact,success);
%    ys = [ys;x(1)];
%    ts = [ts;t];
    if t >= 1 && success% && etanew < param.isopt
        fprintf('Successful approximate solution to %s with p=%3.1e\n',prob.name,p2(1));
        result = 1;
        finaleta = etanew;
        break;
    end
    if itn>=param.itnmax
        fprintf('Maximum iterations reached\n');
        result = 0;
        break;
    end
    
end
g
J
y
% figure(1);
% for i=1:n
%     h(i) = subplot(n+m+1,1,i);
%     plot(ts,xs(:,i));
%     ylabel(strcat('x_',num2str(i)));
%     set(h(i),'ylim',[min(xs(:,i))-sqrt(eps),max(xs(:,i))+sqrt(eps)]);
% end
% for j=1:m
%     h(n+j) = subplot(n+m+1,1,j+n);
%     plot(ts,ys(:,j));
%     ylabel(strcat('y_',num2str(j)));
%     set(h(n+j),'ylim',[min(ys(:,j))-sqrt(eps),max(ys(:,j))+sqrt(eps)]);
% end
% xlabel('t')
% %     h(n+m+1) = subplot(n+m+1,1,m+n+1);
% %     plot(ts,Deltats);
% %     ylabel('\Delta t');
% %     %set(h(n+j),'ylim',[min(ys(:,j))-sqrt(eps),max(ys(:,j))+sqrt(eps)]);
figure(1);
title('dual-degenerate PC-QP');
for i=1:n
    h(i) = subplot(n+m,1,i);
    plot(ts,xs(:,i));
    ylabel(strcat('x(',num2str(i),')'));
    %set(h(i),'ylim',[min(xs(:,i)),max(xs(:,i))]);
    if (min(xs(:,i)) ~= max(xs(:,i)))
        set(h(i),'ylim',[min(xs(:,i)),max(xs(:,i))]);
    end
end
for j=1:m
    h(n+j) = subplot(n+m,1,j+n);
    plot(ts,ys(:,j));
    ylabel(strcat('y(',num2str(j),')'));
    if (min(ys(:,j)) ~= max(ys(:,j)))
        set(h(n+j),'ylim',[min(ys(:,j)),max(ys(:,j))]);
    end
end
xlabel('t');


figure(2);
plot(ts,Etats);
ylabel('\eta');
xlabel('t');
title('dual-degenerate PC-QP');
%linkaxes(h')

%plot(ts,ys,'-.o')
end

function [eta,z] = calceta(x,y,prob,f,g,c,J)
vl = c(prob.me+1:end);

z = g-J'*y;

cmp      = max(0,y(prob.me+1:end).*min(vl,1));
vioL     = min(0,vl);
vioE     = c(1:prob.me);


maxComp  = norm( cmp ,Inf);
maxViol  = norm(    [ vioL; vioE ] ,Inf);
eta = .1*max(maxComp,max(maxViol,norm(z,Inf)));
end

function param = parampath
param = struct('gamma', 0.7, ...
    'gamma1', 0.6, ...
    'itnmax',50000, ...%10000
    'deltat0',.05, ...%.05
    'maxfail',3, ...
    'isopt',0.01,...  %1e-6
    'etamin',10^(-6));
end