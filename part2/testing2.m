% 
% n = randi(5);
% s = randi([2,5]);
% 
% p = randi(12,2,s)-6;
% l = randi(5,n,1);
% 
% mu = 2;
% mu2=1.2;
% 
% lambda = zeros(2*s,1);
% theta = zeros(n*s,1);
% 
% maxAngle=pi/2;


l = [3;2;2];
maxAngle=pi/2;
%p = [-2 0 2 4;4 4 -4 -2];
%p = [-2 -2 -2 -2;4 -3 3 -5];
%p = [0 2 -2;-4 -5 4];
p = [-1 -3 -3 0 3;5 3 -4 5 2];
%p = [5 4 6 4 5;0 2 0.5 -2 -1];
theta = zeros(length(l)*length(p(1,:)),1);
lambda = zeros(length(p(1,:))*2,1);
mu = 2;
mu2 = 1.1;

[theta,eval] = boundConstrainedLagrangian(l,p,theta,lambda,mu2,maxAngle,mu,0.5,0.5)
plotBat(l,theta,p,1);