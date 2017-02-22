function [bool,theta] = checkEnd(l,theta,p)
e = f(l,theta,p);
epsi = 1E-3;
for i = 1:length(l)
    theta(i) = theta(i)+epsi;
    if f(l,theta,p)<e
        bool=0;
        return
    end
    theta(i)=theta(i)-2*epsi;
    if f(l,theta,p)<e
        bool=0;
        return
    end
    theta(i)=theta(i)+epsi;
end
bool = 1;
end
