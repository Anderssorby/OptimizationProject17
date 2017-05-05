function res = gradE(theta,n,s)
%Gradient of the energy function E

res = 2*theta;

for i = 2:s-1
    res(n*(i-1)+1:i*n) = res(n*(i-1)+1:i*n) - theta(n*(i-2)+1:n*(i-1))-theta(n*i+1:n*(i+1));
end

res(1:n) = res(1:n) - theta(n+1:2*n) - theta((s-1)*n+1:end);
res((s-1)*n+1:end) = res((s-1)*n+1:end) - theta(1:n) - theta((s-2)*n+1:(s-1)*n);

end