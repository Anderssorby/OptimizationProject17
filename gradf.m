function res = gradf(l,theta,p)
    n = length(l);
    res = zeros(n,1);
    for i = 1:n
        s1 = 0;
        s2 = 0;
        k1 = 0;
        k2 = 0;
        for j = i:n
            s1 = s1+l(j)*cos(sum(theta(1:j)));
            s2 = s2+l(j)*sin(sum(theta(1:j)));
        end
        for j = 1:n
            k1 = k1+l(j)*cos(sum(theta(1:j)));
            k2 = k2+l(j)*sin(sum(theta(1:j)));
        end
        res(i) = -2*(k1-p(1))*s2+2*(k2-p(2))*s1;
    end
