function alpha = naive_line_search(c1, phi, pk)
phi0 = phi(0);
rho = 0.5;
alpha=1;
% Line search
while phi(alpha) > phi0 + c1*alpha*(pk'*pk)
    alpha = rho*alpha;
end

end
