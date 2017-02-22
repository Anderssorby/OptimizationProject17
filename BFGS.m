function [theta, k] = BFGS(l,p, x0)

max_steps = 1000;

n = length(l);

%B0
Bk = eye(n);
% x0
xk = x0;


gf = gradf(l, xk, p);
k = 0;
c1 = 0.25;
rho = 0.5;
while k < max_steps && norm(gf) > 1E-7
    % This can be more efficient
    pk = Bk\(-gf);
    
    %alpha=1;
    
    % Line search
    %while f(l, xk + alpha*pk, p) > f(l, xk, p) + c1*alpha*(pk'*pk)
    %    alpha = rho*alpha;
    %end
    
    phi = @(alpha) f(l, xk + alpha*pk, p);

    % This could be more efficient
    phiBar = @(alpha) gradf(l, xk + alpha*pk, p)'*pk;
    
    % Line search which satisfy Wolfe conditions
    alpha = line_search(phi, phiBar);
    
    xk = xk + alpha*pk;
    sk = alpha*pk;
    gfkp1 = gradf(l, xk, p);
    yk = gfkp1 - gf;

    Bk = Bk + (yk*yk')/(yk'*sk) - (Bk*sk*sk'*Bk)/(sk'*Bk*sk);

    % finish step
    k = k + 1;
    gf = gfkp1;
end

theta = xk

end
