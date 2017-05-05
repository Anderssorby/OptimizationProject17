
n = randi(9);
s = randi([2,9]);

p = randi(12,2,s)-6;
l = randi(5,n,1);

mu = 0.01;
lambda = zeros(2*s,1);
theta = zeros(n*s,1);

[theta,eval] = equalityConstraintsLagrangian(l,p,theta,lambda,mu,1)
plotArms(l,theta,p,0);