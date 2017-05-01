function [F1,F2] = gradEFFvec(l,theta,n,s,pos)
a = cumsum(theta);
F1 = -cumsum(l.*sin(a),'reverse');
F2 = cumsum(l.*cos(a),'reverse');

F1 = [zeros((pos-1)*n,1);F1;zeros((s-pos)*n,1)];
F2 = [zeros((pos-1)*n,1);F2;zeros((s-pos)*n,1)];

end
