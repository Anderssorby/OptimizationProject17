function res = gradLaGrange(thetavec,lambda,mu,n,s,pmat,l)

res = gradE(thetavec,n,s);

for i = 1:s
    res = res + (mu/2*f(l,thetavec(n*(i-1)+1:n*i),pmat(i,:)) - lambda(i))*gradf(l,thetavec(n*(i-1)+1:n*i),pmat(i,:));
end

end
