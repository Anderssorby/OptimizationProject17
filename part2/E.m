function res = E(thetavec,n,s)
res=0;
for i = 1:s-1
    res = res + norm(thetavec(i*n+1:n*(i+1))-thetavec(n*(i-1)+1:i*n),2)^2;
end
res = res + norm(thetavec(1:n)-thetavec((s-1)*n+1:end),2)^2;
res = res/2;
end
