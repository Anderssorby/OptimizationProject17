function theta = initVals(l,p,theta,maxAngle,lambda,mu)
%Finds initial values for a problem by minimising the augmented lagrangian
%without the energy function E. Uses the gradient box projection method

%l = nx1-dimensional vector with lengths of arm segments
%p = 2xn-dimensional matrix with coordinates for points you want to reach
%theta = s*nx1-dimensional vector with values for all angles
%maxAngle = max allowed angle for the algorithm
%lambda = 2*sx1-dimensional vector with initial values for lambda
%mu = penalty parameter 

%Choose parameters for line search
c1=0.1;
c2=0.9;

%Define s and n
s = length(p(1,:));
n = length(theta)/s;

%convergence tolerance
tolEtha = 0.001;
tolOmega = 0.01;

% start params
omega = 1/mu;
etha = 1/(mu^0.1);

%Projection function
P = @(theta) boxProjection(theta, -maxAngle, maxAngle);

%Gradient descent function with box projection
gd = @(theta,mu,lambda,omega,func,gradFunc) gradientBoxProjection(theta,omega,P,func,gradFunc,c1,c2);

k = 0;
while 1
    
    %Define augmented lagrangian and its gradient
    func = @(theta) laGrangeInit(theta,lambda,mu,n,s,p,l);
    gradFunc = @(theta) gradLaGrangeInit(theta,lambda,mu,n,s,p,l);
    
    % minimise with given tolerances
    theta = gd(theta, mu, lambda,omega,func,gradFunc);

    %Define C-vec
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