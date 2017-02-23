function res = gradf(l,theta,p)
    n = length(l);
    res = zeros(n,1);
    for i = 1:n
        s1 = 0;
        s2 = 0;
        k1 = 0;
        k2 = 0;
        for j = 1:i-1
            a = sum(theta(1:j));
            k1 = k1+l(j)*cos(a);
            k2 = k2+l(j)*sin(a); 
        end
        for j = i:n
            a = sum(theta(1:j));
            k1 = k1+l(j)*cos(a);
            k2 = k2+l(j)*sin(a);
            s1 = s1+l(j)*cos(a);
            s2 = s2+l(j)*sin(a);
        end
        res(i) = -2*(k1-p(1))*s2+2*(k2-p(2))*s1;
end
