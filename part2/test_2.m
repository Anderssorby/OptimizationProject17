

pmat = [5 4 6 4 5; 0 2 0.5 -2 -1];
l = [3; 2; 2];
n = length(l);
s = 5;
theta0 = zeros(n*s,1);
lambda = zeros(2*s,1);
mu2 = 1.2;
mu = 2;
maxAngle = pi/2;

pmat2 = [-1 -3 -3 0 3; 5 3 -4 5 2];
l2 = [];


[theta,eval] = boundConstrainedLagrangian(l,pmat2,theta0,lambda,mu2,maxAngle,mu)

plotBat(l, theta, pmat2)
