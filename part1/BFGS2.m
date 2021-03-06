function [theta,tocvec,fvec] = BFGS2(l,p, x0)
%Start timer
tic

max_steps = 1000;

n = length(l);

%Approximation of the inverse of the hessian
Hk = eye(n);

xk = x0;

%Calculate gradient
gf = gradf(l, xk, p);
k = 0;
c1 = 1E-3;
rho = 0.5;
while 1
    % checking for errors
    if k > max_steps
        error('Surpassed max number of iterations');
    elseif isnan(xk(1))
        error('NaN');
    end
    %Stop condition
    if norm(gf)<1E-7
        [b,xk] = checkEnd(l,xk,p);
        gf = gradf(l, xk, p);
        if b
            break
        end
    end
    % Calculate search direction
    pk = Hk*(-gf);

    phi = @(alpha) f(l, xk + alpha*pk, p);
    phiBar = @(alpha) gradf(l, xk + alpha*pk, p)'*pk;
    
    % Line search which satisfy Wolfe conditions
    alpha = line_search(phi, phiBar);
    
    xk = xk + alpha*pk;
    sk = alpha*pk;
    gfkp1 = gradf(l, xk, p);
    yk = gfkp1 - gf;
    rhok = yk'*sk;
    rhok=1/rhok;

    if rhok <= 0
        % This will create a not positive definite matrix
        % and ruin our method. Therefore do not update the matrix.
    else
        %Update approximation matrix
        Hk = (eye(n)-rhok*sk*yk')*Hk*(eye(n)-rhok*yk*sk')+rhok*sk*sk';
    end
    
    % finish step
    k = k + 1;
    gf = gfkp1;
    %Update error vector
    fvec(k) = sqrt(f(l,xk,p));
    %Update timer vecor
    tocvec(k) = toc;
    
end

theta = xk;

end
