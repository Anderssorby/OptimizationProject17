function bool = reachAPoint(l,p,theta,maxAngle,mu)
% WORK IN PROGRESS
%x0
%lambda0
bool=0;
%reachTol = 1E-4;
%ctol = 0.0001;

% if ~checkReachable(l,pmat,reachTol)
%     error('Not reachable');
% end

n = length(theta);

%convergence tolerance

tolOmega = 0.001;
tolEtha = 0.001;

% start params
%mu = 10;
omega = 1/mu;
etha = 1/mu^0.1;

P = @(theta) boxProjection(theta, -maxAngle, maxAngle);

%gradBCL = @(theta, lambda, mu, cvec, gradCvec) gradf(l, p, theta)-lambda*cvec'+mu*sum(cvec.*gradCvec);
func = @(theta) f(l,theta,p);
gradFunc = @(theta) gradf(l,theta,p);
gd = @(theta,omega) BFGSbox(theta,omega,P,func,gradFunc);
%gf = gradLaGrange(thetak,lambda,mu,n,s,p,l);
k = 0;
while 1
    % solve such that modified KKT holds for less than tolerance
    theta = gd(theta,omega)
    if norm(bigEff(l,theta,n)-p)<tolEtha
        if norm(bigEff(l,theta,n)-p)<tolEtha && norm(theta-P(theta-gradFunc(theta))) <= tolOmega
            bool=1;
            % Acceptable solution found
            break;
        end
        % Update multipliers tighten tolerances
        omega = omega/mu;
        etha = etha/mu^0.9;
    else
        mu = 100*mu;
        etha = 1/mu^0.1;
        omega = 1/mu;
    end
    k = k + 1;
end

end