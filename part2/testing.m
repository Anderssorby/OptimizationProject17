
n = randi(9);
s = randi([2,9]);

p = randi(12,2,s)-6;
l = randi(5,n,1);

mu = 1;
lambda = zeros(2*s,1);
theta = zeros(n*s,1);

[theta,eval] = augLagr(l,p,theta,lambda,mu,1)
plotBat(l,theta,p);