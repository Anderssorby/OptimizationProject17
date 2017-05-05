function L = laGrange(thetavec,lambda,mu,n,s,p,l)
%Calculates the augmented Lagrangian with a given theta, lambda, mu and p
L = E(thetavec,n,s);

for i = 1:s
    L = L - lambda(2*i-1:2*i)'*(bigEff(l,thetavec(n*(i-1)+1:n*i),n)-p(:,i)) + mu*(f(l,thetavec(n*(i-1)+1:n*i),p(:,i)));
end

end
    