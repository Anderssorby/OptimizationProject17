function L = laGrangeInit(thetavec,lambda,mu,n,s,pmat,l)
%Augmented Lagrangian without energy funciton E

L = 0;

for i = 1:s
    L = L - lambda(2*i-1:2*i)'*(bigEff(l,thetavec(n*(i-1)+1:n*i),n)-pmat(:,i)) + mu*(f(l,thetavec(n*(i-1)+1:n*i),pmat(:,i)));
end

end
    