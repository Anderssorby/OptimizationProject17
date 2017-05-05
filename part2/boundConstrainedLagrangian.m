function [theta,eval] = boundConstrainedLagrangian(l,pmat,theta,lambda,mu2,maxAngle,mu,tolOmega,tolEtha)

%Find initial values for theta
theta = initVals(l,pmat,theta,maxAngle,lambda,mu2); 

%Choose parameters for line search
c1=0.1;
c2=0.9;

%Define s and n
s = length(pmat(1,:));
n = length(theta)/s;


% start params
omega = 1/mu;
etha = 1/(mu^0.1);

%Projection function
P = @(theta) boxProjection(theta, -maxAngle, maxAngle);

%line search algorithm
gd = @(theta,mu,lambda,omega,func,gradFunc) gradientBoxProjection(theta,omega,P,func,gradFunc,c1,c2);

k = 0;
while 1
    %Define augmented lagrangian and its gradient
    func = @(theta) laGrange(theta,lambda,mu,n,s,pmat,l);
    gradFunc = @(theta) gradLaGrange(theta,lambda,mu,n,s,pmat,l);
    
    % minimise with given tolerances
    theta = gd(theta, mu, lambda,omega,func,gradFunc);
    
    %Make C-vector
    cvec = zeros(s,1);
    for i = 1:s
        cvec(2*i-1:2*i) = bigEff(l,theta((i-1)*n+1:i*n),n)-pmat(:,i);
    end

    if norm(cvec) <= etha
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
eval = E(theta,n,s); %Calculate E(theta)
end