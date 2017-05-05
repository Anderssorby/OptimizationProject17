function [theta,eval] = equalityConstraintsLagrangian(l,pmat,theta,lambda,mu,tol)

%Tolerance for C-vec
ctol = 1E-5;

%Choose parameters for line search
c1=0.1;
c2=0.9;

if ~checkAnnulus(l,pmat)
    error('Not reachable');
end

s = length(pmat(1,:));
n = length(theta)/s;

%Slacking/tightening factors for mu and tau
ksi = 2;
gamma = 0.8;

k=1;

while 1
    func = @(theta) laGrange(theta,lambda,mu,n,s,pmat,l);
    gradFunc = @(theta) gradLaGrange(theta,lambda,mu,n,s,pmat,l);
    theta = BFGS(theta,tol,func,gradFunc,c1,c2)
    %theta = gradientDescent(theta,tol,func,gradFunc,c1)
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