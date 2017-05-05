function theta = initVals(l,p,theta,maxAngle,lambda,mu)

c1=0.1;
c2=0.9;

s = length(p(1,:));
n = length(theta)/s;

%convergence tolerance
tolEtha = 0.001;
tolOmega = 0.01;

% start params
%mu = 10;
omega = 1/mu;
etha = 1/(mu^0.1);

P = @(theta) boxProjection(theta, -maxAngle, maxAngle);

%gradBCL = @(theta, lambda, mu, cvec, gradCvec) gradE(l, pmat, theta)-lambda*cvec'+mu*sum(cvec.*gradCvec);

gd = @(theta,mu,lambda,omega,func,gradFunc) gradientBoxProjection(theta,omega,P,func,gradFunc,c1,c2);
%gf = gradLaGrange(thetak,lambda,mu,n,s,p,l);
k = 0;
while 1
    func = @(theta) laGrangeInit(theta,lambda,mu,n,s,p,l);
    gradFunc = @(theta) gradLaGrangeInit(theta,lambda,mu,n,s,p,l);
    % solve such that modified KKT holds for less than tolerance
    %try
    theta = gd(theta, mu, lambda,omega,func,gradFunc)
    %catch ME
    %    error('Could not hit all points with initVals')
    %end
    %if checkConstraints(l,theta,pmat,n,s,ctol)
    %    break;
    %end
    cvec = zeros(s,1);
    for i = 1:s
        cvec(2*i-1:2*i) = bigEff(l,theta((i-1)*n+1:i*n),n)-p(:,i);
    end

    if norm(cvec) <= etha
        % test convergence

        if norm(cvec) <= tolEtha && norm(theta-P(theta-gradLaGrangeInit(theta,lambda,mu,n,s,p,l))) <= tolOmega
            % Acceptable solution found
            break;
        end
        % Update multipliers tighten tolerances
        lambda = lambda - mu*cvec;
        %mu = mu;
        etha = etha/(mu^0.9);
        omega = omega/mu;
    else
        % Increase penalty parameter, tighten tolerances
        mu = 2*mu;
        etha = 1/(mu^0.1);
        omega = 1/mu;
    end

    k = k + 1;
end

end