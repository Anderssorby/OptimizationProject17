function res = f(l,theta,p)
%l = Vector with lengths of joints
%theta = Vector with angles between joints
%p = point we want the distance to
    l=l';
    a = cumsum(theta);
    s1 = l*cos(a);
    s2 = l*sin(a);
    res = (s1-p(1))^2+(s2-p(2))^2; %Calculate the squared 2-norm of F-p
end