function alpha = zoom(phi, phiBar, constants, a_lo, a_hi)
% Assuming that the interval [a_lo, a_hi] contains step length
% that satisfies the Wolfe conditions
% That means that a_lo satisifes while a_hi doesn't satisfy the Armijio rule
% alpha = Step length for next step of our algorithm
% phi = function phi(a) = f(x+a*p)
% phiBar = the derivative of phi(a) wrt a, multiplied with p_k


c1 = constants(1);
c2 = constants(2);

alpha = a_lo;
prevalpha = alpha;

phiBar0 = phiBar(0);
phi0 = phi(0);

c = 0;

while c < 10

    % Interpolate to find a test step length alpha

    % Minimizing the interpolation on [a_lo, a_hi]
    % This is our trial step length
    alpha = interpolateQuad(a_lo, a_hi, phi(a_lo), phiBar(a_lo), phi(a_hi));
   
    % Naive interpolation
    if alpha <= 0
        alpha = (a_hi + a_lo)/2;
    end

    phia = phi(alpha);
    
    if phia > phi0 + c1*alpha*phiBar0 || phia >= phi(a_lo)
        % Alpha does also not satisfy the Armijio rule
        % or a_lo is a better alpha than alpha
        % Our alpha is an improvement of a_hi
        a_hi = alpha;
    else
        % alpha satisfies the Armijo rule and is an improvement to a_lo
        phiBara = phiBar(alpha);
        if abs(phiBara) <= -c2*phiBar0
           % aplpha satisfies the second Wolfe condition
           % and we are done
           break
        end
        if phiBara*(a_hi - a_lo) >= 0
           % Change the interval
           a_hi = a_lo;
        end
        a_lo = alpha;
    end
    c = c + 1;

    h = a_hi-a_lo;
    if h < 1E-4
        % Don't let h be ridiculously small
        alpha = a_lo;
        break
    end


    % Look out for errors
    if alpha <= 0
        error('Zoom: Alpha less than 0')
    end
    prevalpha = alpha;
end


end