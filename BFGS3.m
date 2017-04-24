function [theta,tocvec,fvec] = BFGS2(l,p, x0)
tic
max_steps = 1000;

n = length(l);

%B0
Hk = eye(n);
% x0
xk = x0;


gf = gradf(l, xk, p);
k = 0;
c1 = 1E-3;
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
    pk = Hk*(-gf);
    phi = @(alpha) f(l, xk + alpha*pk, p);

    % This could be more efficient
    phiBar = @(alpha) gradf(l, xk + alpha*pk, p)'*pk;
    
    % Line search which satisfy Wolfe conditions
    alpha = line_search(phi, phiBar);
    xk = xk + alpha*pk;
    sk = alpha*pk;
    gfkp1 = gradf(l, xk, p);
    yk = gfkp1 - gf;
    rhok = yk'*sk;
    rhok=1/rhok;
    %Bk = Bk + (yk*yk')/(yk'*sk) - (Bk*sk*sk'*Bk)/(sk'*Bk*sk);
    Hk = (eye(n)-rhok*sk*yk')*Hk*(eye(n)-rhok*yk*sk')+rhok*sk*sk';
    
    % finish step
    k = k + 1;
    gf = gfkp1;
    fvec(k) = sqrt(f(l,xk,p));
    tocvec(k) = toc;
    
end

theta = xk;
%disp(k)
end