function res = gradLaGrange(thetavec,lambda,mu,n,s,pmat,l)

res = gradE(thetavec,n,s);

% for i = 1:s
%     for j = 1:2
%         res = res - 
%     end
%     res = res + (mu/2*f(l,thetavec(n*(i-1)+1:n*i),pmat(:,i)) - lambda(i))*gradfvec(l,thetavec(n*(i-1)+1:n*i),pmat(:,i),i,s);
% end

for i = 1:s
    [FF1,FF2] = gradEFFvec(l,thetavec(n*(i-1)+1:n*i),n,s,i);
    res = res + mu*gradfvec(l,thetavec(n*(i-1)+1:n*i),pmat(:,i),i,s) - lambda(2*i-1)*FF1-lambda(2*i)*FF2;
end

end
