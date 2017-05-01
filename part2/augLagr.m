function [theta,eval] = augLagr(l,pmat,theta,lambda,mu,tol)
reachTol = 1E-4;
ctol = 1E-5;

if ~checkReachable(l,pmat,reachTol)
    error('Not reachable');
end

s = length(pmat(1,:));
n = length(theta)/s;

ksi = 3;
gamma = 0.8;

gd = @(theta,tol,mu,lambda) BAFGS(l,pmat,theta,tol,mu,lambda,n,s);
k=1;

while 1
    theta = gd(theta,tol,mu,lambda)
    if checkConstraints(l,theta,pmat,n,s,ctol)
        break;
    end
    cvec = zeros(2*s,1);
    for i = 1:s
        cvec(2*i-1:2*i) = bigEff(l,theta(n*(i-1)+1:n*i),n)-pmat(:,i);
        %cvec(i) = f(l,theta((i-1)*n+1:i*n),pmat(:,i));
    end
    lambda = lambda - mu*cvec;
    mu = ksi*mu;
    tol = gamma*tol;
    k=k+1;
    if k >= 1000
        error('Max number of iterations surpassed (augLagr)');
    end
end
eval = E(theta,n,s);
end