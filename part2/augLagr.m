function theta = augLagr(l,pmat,theta,lambda,mu,tol)
reachTol = 1E-4;
ctol = 1;

if ~checkReachable(l,pmat,reachTol)
    error('Not reachable');
end

s = length(pmat(1,:));
n = length(theta)/s;

ksi = 3;
gamma = 0.8;

gd = @(theta,tol,mu,lambda) gradientDescentAugment(l,pmat,theta,tol,mu,lambda,n,s);
k=1;

while 1
    theta = gd(theta,tol,mu,lambda)
    if checkConstraints(l,theta,pmat,n,s,ctol)
        break;
    end
    cvec = zeros(s,1);
    for i = 1:s
        cvec(i) = f(l,theta((i-1)*n+1:i*n),pmat(:,i));
    end
    lambda = lambda - mu*cvec;
    mu = ksi*mu;
    tol = gamma*tol;
    k=k+1;
    if k >= 1000
        error('Max number of iterations surpassed (augLagr)');
    end
end
end