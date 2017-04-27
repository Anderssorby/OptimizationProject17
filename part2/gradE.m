function res = gradE(thetavec,n,s)

res = 2*thetavec;

for i = 2:s-1
    res(n*(i-1)+1:i*n) = res(n*(i-1)+1:i*n) - thetavec(n*(i-2)+1:n*(i-1))-thetavec(n*i+1:n*(i+1));
end

res(1:n) = res(1:n) - thetavec(n+1:2*n) - thetavec((s-1)*n+1:end);
res((s-1)*n+1:end) = res((s-1)*n+1:end) - thetavec(1:n) - thetavec((s-2)*n+1:(s-1)*n);

end