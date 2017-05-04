function [theta,eval] = boundConstrainedLagrangian(l,pmat,theta,lambda,mu2,maxAngle,mu)
% WORK IN PROGRESS
%x0
%lambda0
theta = initVals(l,pmat,theta,maxAngle,lambda,mu2);


c1=0.1;
c2=0.9;
% if ~checkReachable(l,pmat,reachTol)
%     error('Not reachable');
% end

s = length(pmat(1,:));
n = length(theta)/s;

%convergence tolerance
tolEtha = 0.4;
tolOmega = 0.5;

% start params
%mu = 10;
omega = 1/mu;
etha = 1/(mu^0.1);

P = @(theta) boxProjection(theta, -maxAngle, maxAngle);

%gradBCL = @(theta, lambda, mu, cvec, gradCvec) gradE(l, pmat, theta)-lambda*cvec'+mu*sum(cvec.*gradCvec);

gd = @(theta,mu,lambda,omega,func,gradFunc) BFGSbox(theta,omega,P,func,gradFunc,c1,c2);
%gf = gradLaGrange(thetak,lambda,mu,n,s,p,l);
k = 0;
while 1
    func = @(theta) laGrange(theta,lambda,mu,n,s,pmat,l);
    gradFunc = @(theta) gradLaGrange(theta,lambda,mu,n,s,pmat,l);
    % solve such that modified KKT holds for less than tolerance
    theta = gd(theta, mu, lambda,omega,func,gradFunc)
    %if checkConstraints(l,theta,pmat,n,s,ctol)
    %    break;
    %end
    cvec = zeros(s,1);
    for i = 1:s
        cvec(2*i-1:2*i) = bigEff(l,theta((i-1)*n+1:i*n),n)-pmat(:,i);
    end

    if norm(cvec) <= etha
        disp('ja')
        disp(lambda)
        disp(cvec)
        % test convergence

        if norm(cvec) <= tolEtha && norm(gradFunc(theta)) <= tolOmega
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
eval = E(theta,n,s);
end