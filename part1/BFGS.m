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
while 1

    if k > max_steps
        error('Surpassed max number of iterations');
    elseif isnan(xk(1))
        error('NaN');
    end
    if norm(gf)<1E-7
        [b,xk] = checkEnd(l,xk,p);
        gf = gradf(l, xk, p);
        if b
            break
        end
    end
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

    yksk = yk'*sk;
    if yksk <= 0
        % This will create a not positive definite matrix
        % and ruin our method. Therefore do not update the matrix.
    else
        Bk = Bk + (yk*yk')/yksk - (Bk*sk*sk'*Bk)/(sk'*Bk*sk);
    end

    % finish step
    k = k + 1;
    gf = gfkp1;
end

theta = xk

end
