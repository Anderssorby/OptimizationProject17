function res = bigEff(l,theta,n)
    s1 = 0;
    s2 = 0;
    for i = 1:min(n,length(l))
        s1 = s1+l(i)*cos(sum(theta(1:i)));
        s2 = s2+l(i)*sin(sum(theta(1:i)));
    end
    res=[s1,s2];
end