function [alpha] = line_search(phi, phiBar,c1,c2)


% how much to scale alpha in each step
inc = 1.5;

% Store values
prevphi = 0;
prevalpha = 0;
alpha = 1;


phi0 = phi(0);
phiBar0 = phiBar(0);

constants = [c1, c2];

i = 0;
while i < 10

    if alpha <= 0
        alpha
        error('Alpha less than 0')
    end

    phia = phi(alpha);
    if phia > phi0 + c1*alpha*phiBar0 || (phia > prevphi && i > 1)
        % alpha does not satisfy Armijo rule and we can use 
        % it as an upper bound in zoom
        alpha = zoom(phi, phiBar, constants, prevalpha, alpha);
        break
    end
    
    phiBara = phiBar(alpha);
    
    if abs(phiBara) <= -c2*phiBar0
        % Second Wolfe condition satisfied
        break
    end

    if phiBara >= 0
        % We are overstepping
        alpha = zoom(phi, phiBar, constants, alpha, prevalpha);
        break
    end
    
    prevalpha = alpha;
    % choose alpha
    alpha = alpha*inc;


    i = i + 1;
end

end