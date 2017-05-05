function res = bigEff(l,theta,n)
%Calculates position of n'th joint of the robot arm
l=l';
a = cumsum(theta);
m = min(n,length(l));
l = l(1:m);
a = a(1:m);
res=[l*cos(a);l*sin(a)];
end