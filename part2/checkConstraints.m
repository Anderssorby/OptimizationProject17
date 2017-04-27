function bool = checkConstraints(l,theta,pmat,n,s,ctol)
bool = false;
su = 0;
for i = 1:s
    su = su + f(l,theta(n*(i-1)+1:n*i),pmat(:,i));
end
if su < ctol
    bool = true;
end
end