function res = f(l,theta,p)
%l = Vector with lengths of joints
%theta = Vector with angles between joints
%p = point we want the distance to
    s1 = 0;
    s2 = 0;
    for i = 1:length(l)
        a = sum(theta(1:i)); %sum of angles up to index i
        s1 = s1+l(i)*cos(a); %1st coordinate of the vector F
        s2 = s2+l(i)*sin(a); %2nd coordinate of the vector F
    end
    res = (s1-p(1))^2+(s2-p(2))^2; %Calculate the squared 2-norm of F-p
end
