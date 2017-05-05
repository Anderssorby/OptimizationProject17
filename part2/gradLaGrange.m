function res = gradLaGrange(thetavec,lambda,mu,n,s,pmat,l)
%Gradient of the augmented Lagrangian

res = gradE(thetavec,n,s);



for i = 1:s
    [c1,c2] = gradCvec(l,thetavec(n*(i-1)+1:n*i),n,s,i);
    res = res + mu*gradfvec(l,thetavec(n*(i-1)+1:n*i),pmat(:,i),i,s) - lambda(2*i-1)*c1-lambda(2*i)*c2;
end

end
