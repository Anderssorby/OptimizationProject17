function reachable = checkReachable(l,pmat,tol)

reachable = 1;
theta0 = rand(length(l),1);

for i = 1:length(pmat(1,:))
    [t,theta,fvec]=BFGS2(l,pmat(:,i),theta0);
    if fvec(end)>tol
        reachable=0;
        return;
    end
end

end