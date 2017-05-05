function theta = gradientBoxProjection(theta0,omega,P,func,gradFunc,c1,c2)

max_steps = 1000;

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
    if norm(gf)<omega
        break;
    end
    
    
    % Calculate search direction
    pk = -gf;

    %phi = @(alpha) laGrange(thetak + alpha*pk,lambda,mu,n,s,p,l);
    phi = @(alpha) func(thetak + alpha*pk);
    phiBar = @(alpha) gradFunc(thetak + alpha*pk)'*pk;
    
    % Line search which satisfy Wolfe conditions
    alpha = line_search(phi, phiBar,c1,c2);
    
    thetak = P(thetak + alpha*pk);
    gf = gradFunc( thetak);

    
    % finish step
    k = k + 1;

    
end

theta = thetak;

end