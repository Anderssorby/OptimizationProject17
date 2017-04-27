function bool = checkConstraints(theta,pmat,s,ctol)
bool = false;
su = 0;
for i = 1:s
    su = su + f(l,theta(s*(i-1)+1:s*i),pmat(:,i));
end
if su < ctol
    bool = true;
end
end