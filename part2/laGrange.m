function L = laGrange(thetavec,lambda,mu,n,s,pmat,l)

L = E(thetavec,n,s);

for i = 1:s
    L = L - lambda(i)*f(l,thetavec(n*(i-1)+1:n*i),pmat(i,:)) + mu/2*(f(l,thetavec(n*(i-1)+1:n*i),pmat(i,:)))^2;
end

end
    