function theta = boundConstrainedLagrangian(l,pmat,theta,lambda,mu,tol, maxAngle)
% WORK IN PROGRESS
%x0
%lambda0

reachTol = 1E-4;
ctol = 1;

if ~checkReachable(l,pmat,reachTol)
    error('Not reachable');
end

s = length(pmat(1,:));
n = length(theta)/s;

%convergence tolerance
tolEtha = 0.1;
tolOmega = 0.1;

% start params
mu = 10;
omega = 1/mu;
etha = 1/(mu^0.1);

P = @(theta) boxProjection(theta, -maxAngle, maxAngle);

gradBCL = @(theta, lambda, mu, cvec, gradCvec) gradE(l, pmat, theta)-lambda*cvec'+mu*sum(cvec.*gradCvec);

gd = @(theta,tol,mu,lambda) gradientDescentAugment(l,pmat,theta,tol,mu,lambda,n,s);

k = 0;
while 1
    % solve such that modified KKT holds for less than tolerance
    theta = gd(theta, omega, mu, lambda);
    %if checkConstraints(l,theta,pmat,n,s,ctol)
    %    break;
    %end
    cvec = zeros(s,1);
    for i = 1:s
        cvec(i) = f(l,theta((i-1)*n+1:i*n),pmat(:,i));
    end

    if norm(cvec) <= etha
        % test convergence

        if norm(cvec) <= tolEtha && norm(theta-P(theta-gradBCL())) <= tolOmega
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
        lambda = lambda;
        mu = 100*mu;
        etha = 1/(mu^0.1);
        omega = 1/mu;
    end

    k = k + 1;
end
end
