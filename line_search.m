function [alpha] = line_search(phi, phiBar)
% Algorithm 3.5
amax = 1;
% alpha0
c1 = 1E-4;
c2 = 0.1;

% how much to increment alpha in each step
inc = 0.01;

prevphi = 0;

prevalpha = 0;
alpha = 0.001;


phi0 = phi(0);
phiBar0 = phiBar(0);
i = 1;
while i < 1000

    phia = phi(alpha);
    if phia > phi0 + c1*alpha*phiBar0 || (phia > prevphi && i > 1)
        alpha = zoom(prevalpha, alpha);
        break
    end
    
    phiBara = phiBar(alpha);
    
    if abs(phiBara) <= -c2*phiBar0
        break
    end

    if phiBara >= 0
        alpha = zoom(alpha, prevalpha);
        break
    end

    % choose alpha
    alpha = alpha + inc;


    i = i + 1;
end
end
