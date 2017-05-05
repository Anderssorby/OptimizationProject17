function res = E(theta,n,s)
%The energy function E
res=0;
for i = 1:s-1
    res = res + norm(theta(i*n+1:n*(i+1))-theta(n*(i-1)+1:i*n),2)^2;
end
res = res + norm(theta(1:n)-theta((s-1)*n+1:end),2)^2;
res = res/2;
end
