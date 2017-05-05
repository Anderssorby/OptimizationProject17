function [theta,eval] = boundConstrainedLagrangian(l,p,theta,lambda,mu2,maxAngle,mu,tolOmega,tolEta)
%l = nx1-dimensional vector with lengths of arm segments
%p = 2xn-dimensional matrix with coordinates for points you want to reach
%theta = s*nx1-dimensional vector with values for all angles
%lambda = 2*sx1-dimensional vector with initial values for lambda
%mu2 = penalty parameter for finding initial values
%maxAngle = max allowed angle for the algorithm
%mu = penalty parameter the gradient projection method
%tolOmega = tolerance for the norm of the gradient in the gradient projection method
%tolEta = tolerance for the norm of the C-vec in the gradient projection
%method

%eval = Value of E(theta) for final value of theta




%Find initial values for theta
theta = initVals(l,p,theta,maxAngle,lambda,mu2); 

%Choose parameters for line search
c1=0.1;
c2=0.9;

%Define s and n
s = length(p(1,:));
n = length(theta)/s;


% start params
omega = 1/mu;
eta = 1/(mu^0.1);

%Projection function
P = @(theta) boxProjection(theta, -maxAngle, maxAngle);

%Gradient descent function with box projection
gd = @(theta,mu,lambda,omega,func,gradFunc) gradientBoxProjection(theta,omega,P,func,gradFunc,c1,c2);

k = 0;
while 1
    %Define augmented lagrangian and its gradient
    func = @(theta) laGrange(theta,lambda,mu,n,s,p,l);
    gradFunc = @(theta) gradLaGrange(theta,lambda,mu,n,s,p,l);
    
    % minimise with given tolerances
    theta = gd(theta, mu, lambda,omega,func,gradFunc);
    
    %Make C-vector
    cvec = zeros(s,1);
    for i = 1:s
        cvec(2*i-1:2*i) = bigEff(l,theta((i-1)*n+1:i*n),n)-p(:,i);
    end

    if norm(cvec) <= eta
        % test convergence

        if norm(cvec) <= tolEta && norm(gradFunc(theta)) <= tolOmega
            % Acceptable solution found
            break;
        end
        % Update multipliers tighten tolerances
        lambda = lambda - mu*cvec;
        %mu = mu;
        eta = eta/(mu^0.9);
        omega = omega/mu;
    else
        % Increase penalty parameter, tighten tolerances
        mu = 2*mu;
        eta = 1/(mu^0.1);
        omega = 1/mu;
    end

    k = k + 1;
end
eval = E(theta,n,s); %Calculate E(theta)
end