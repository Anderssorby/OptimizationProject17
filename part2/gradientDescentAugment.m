function [theta] = gradientDescentAugment(l,pmat,theta0,tol,mu,lambda,n,s)
% % %theta = vector of angles between the different joints
% % %tocvec = Vector with the time used until each iteration of the method
% % %fvec = Vector with distance from the point p for each iteration
% % 
% % %l = Vector with lengths for the joints of the robot
% % %p = Coordinates for the point we are trying to reach
% % %x0 = Initial values for theta

%tic %start timing 
theta = theta0; %Initiate theta
count = 0; %Initiate counter
c1 = 0.0001; %Choose c1
gl = gradLaGrange(theta,lambda,mu,n,s,pmat,l); %Find gradient at initial position
rho = 0.5; %Define slacking factor
while true
    count = count + 1;
    
    %Backtracking method
    a=1; %Start with alpha = 1
    while laGrange(theta-a*gl,lambda,mu,n,s,pmat,l)>laGrange(theta,lambda,mu,n,s,pmat,l)-c1*a*(gl'*gl)
        a = rho*a; %Decrease alpha when not satisfying Armijo conditions
    end
    
    theta = theta - gl*a; %Find next step for theta
    gl=gradLaGrange(theta,lambda,mu,n,s,pmat,l); %Calculate gradient in new position
    %fvec(count) = sqrt(f(l,theta,p)); %Find position from this step to p
    %tocvec(count) = toc; %Save time until this iteration
    
    %Check for end
    if count > 1000
        error('Surpassed max number of iterations'); %Quit if too many steps have been used
    end
    if norm(gl)<tol
        break
%         [b,theta] = checkEnd(l,theta,p); %Is the endpoint a maximum or a sadle point?
%         if b
%             break %Break if you have reached a minimum
%         end
    end
end
end