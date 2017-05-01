function theta = BFGSbox(l,p,theta0,mu,lambda,n,s,maxAngle,omega)

max_steps = 5000;

P = @(theta) boxProjection(theta, -maxAngle, maxAngle);

%Approximation of the inverse of the hessian
Hk = eye(n*s);

thetak = theta0;

%Calculate gradient
gf = gradLaGrange(thetak,lambda,mu,n,s,p,l);
k = 0;
c1 = 1E-3;
rho = 0.5;
while 1
    % checking for errors
    if k > max_steps
        error('Surpassed max number of iterations');
    elseif isnan(thetak(1))
        error('NaN');
    end
    %Stop condition
    if norm(thetak-P(thetak-gf))<omega
        %[b,thetak] = checkEnd(l,thetak,p);
        %gf = gradLaGrange(l, thetak, p);
        %if b
            break
        %end
    end
    % Calculate search direction
    pk = Hk*(-gf);

    phi = @(alpha) laGrange(thetak + alpha*pk,lambda,mu,n,s,p,l);
    phiBar = @(alpha) gradLaGrange(thetak + alpha*pk,lambda,mu,n,s,p,l)'*pk;
    
    % Line search which satisfy Wolfe conditions
    alpha = line_search(phi, phiBar);
    
    thetak = thetak + alpha*pk;
    sk = alpha*pk;
    gfkp1 = gradLaGrange( thetak,lambda,mu,n,s, p,l);
    yk = gfkp1 - gf;
    rhok = yk'*sk;
    rhok=1/rhok;

    if rhok <= 0
        % This will create a not positive definite matrix
        % and ruin our method. Therefore do not update the matrix.
    else
        %Update approximation matrix
        Hk = (eye(n*s)-rhok*sk*yk')*Hk*(eye(n*s)-rhok*yk*sk')+rhok*sk*sk';
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