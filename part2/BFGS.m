function theta = BFGS(theta0,tol,func,gradFunc,c1,c2)

max_steps = 1000;
n = length(theta0);
%Approximation of the inverse of the hessian
Hk = eye(n);

thetak = theta0;

%Calculate gradient
gf = gradFunc(thetak);
k = 0;

while 1
    % checking for errors
    if k > max_steps
        error('Surpassed max number of iterations');
    elseif isnan(thetak(1))
        error('NaN');
    end
    %Stop condition
    if norm(gf)<tol
        %[b,thetak] = checkEnd(l,thetak,p);
        %gf = gradLaGrange(l, thetak, p);
        %if b
            break
        %end
    end
    % Calculate search direction
    pk = Hk*(-gf);

    phi = @(alpha) func(thetak + alpha*pk);
    phiBar = @(alpha) gradFunc(thetak + alpha*pk)'*pk;
    
    % Line search which satisfy Wolfe conditions
    alpha = line_search(phi, phiBar,c1,c2);
    
    thetak = thetak + alpha*pk;
    sk = alpha*pk;
    gfkp1 = gradFunc( thetak);
    yk = gfkp1 - gf;
    rhok = yk'*sk;
    rhok=1/rhok;

    if rhok <= 0
        % This will create a not positive definite matrix
        % and ruin our method. Therefore do not update the matrix.
    else
        %Update approximation matrix
        Hk = (eye(n)-rhok*sk*yk')*Hk*(eye(n)-rhok*yk*sk')+rhok*(sk*sk');
    end
    
    % finish step
    k = k + 1;
    gf = gfkp1;
    %Update error vector
    %fvec(k) = sqrt(laGrange(thetak,lambda,mu,n,s,p,l));
    %Update timer vecor
    %tocvec(k) = toc;
    
end

theta = thetak;

end