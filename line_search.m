function [alpha] = line_search(phi, phiBar)
% Algorithm 3.5 

% alpha0
amax = 1;

c1 = 1E-4;
c2 = 0.1;

% how much to scale alpha in each step
inc = 1.01;

% Store values
prevphi = 0;
prevalpha = 0;
alpha = 0.001;


phi0 = phi(0);
phiBar0 = phiBar(0);

constants = [c1, c2];

i = 1;
while i < 1000

    if alpha <= 0
        alpha
        error('Alpha less than 0')
    end

    phia = phi(alpha);
    if phia > phi0 + c1*alpha*phiBar0 || (phia > prevphi && i > 1)
        alpha = zoom(phi, phiBar, constants, prevalpha, alpha);
        break
    end
    
    phiBara = phiBar(alpha);
    
    if abs(phiBara) <= -c2*phiBar0
        break
    end

    if phiBara >= 0
        alpha = zoom(phi, phiBar, constants, alpha, prevalpha);
        break
    end
    
    prevalpha = alpha;
    % choose alpha
    alpha = alpha*inc;


    i = i + 1;
end
end
