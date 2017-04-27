function res = gradf(l,theta,p)
%l: the arm lengths
%theta: the current angles
%p: point to aproximate
    a=cumsum(theta);
    s1=cumsum(l.*cos(a),'reverse');
    s2=cumsum(l.*sin(a),'reverse');
    res=-(s1(1)-p(1))*s2+(s2(1)-p(2))*s1;
end
