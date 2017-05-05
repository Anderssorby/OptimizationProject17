function [c1,c2] = gradCvec(l,theta,n,s,pos)
%Returns the gradient of two equality constraints. We need pos, n and s so
%we can place it in the gradient vector where most indices will be zero

a = cumsum(theta);
c1 = -cumsum(l.*sin(a),'reverse');
c2 = cumsum(l.*cos(a),'reverse');

c1 = [zeros((pos-1)*n,1);c1;zeros((s-pos)*n,1)];
c2 = [zeros((pos-1)*n,1);c2;zeros((s-pos)*n,1)];

end
