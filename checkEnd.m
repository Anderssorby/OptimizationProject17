function [bool,theta] = checkEnd(l,theta,p) %Checks if end is a min point
%bool = 1 if endpoint is a min-point
%theta = angles between joints
%l = vector of lengths of joints
%p = point we try to reach

e = f(l,theta,p); % value of f at endpoint
epsi = 1E-3; %Choose some epsilon>0

%Test values of f for close values of theta
for i = 1:length(l)
    theta(i) = theta(i)+epsi;
    if f(l,theta,p)<e %This means we are not in a min point
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
bool = 1; %All close values were larger than e. We are in a min point
end
