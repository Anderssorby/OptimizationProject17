function res = gradfvec(l,theta,p,pos,s)
%l: the arm lengths
%theta: the current angles
%p: point to aproximate
    a=cumsum(theta);
    s1=cumsum(l.*cos(a),'reverse');
    s2=cumsum(l.*sin(a),'reverse');
    res=-(s1(1)-p(1))*s2+(s2(1)-p(2))*s1;
    n = length(res);
    res = [zeros((pos-1)*n,1);res;zeros((s-pos)*n,1)];
end
