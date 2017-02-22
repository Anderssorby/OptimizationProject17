function alpha = zoom(phi, phiBar, constants, a_lo, a_hi)
% Assuming that the interval [a_lo, a_hi] contains step length
% that satisfies the strong Wolfe conditions

c1 = constants(1);
c2 = constants(2);

alpha = a_lo;

phiBar0 = phiBar(0);
phi0 = phi(0);

c = 0;

while c < 1000

    % This shouldn't happen
    if a_lo > a_hi
        error('Zoom: a_lo > a_hi')
    end

    % Interpolate to find a test step length alpha

    % Minimizing the interpolation on [a_lo, a_hi]
    % This is our trial step length
    %alpha = (phiBar(a_lo)*a_hi^2)/(2*(phi(a_hi) - phi(a_lo)- phiBar(a_lo)*a_hi));
   
    % Naive interpolation
    alpha = (a_hi + a_lo)/2;

    a_lo
    

    phia = phi(alpha);
    
    if phia > phi0 + c1*alpha*phiBar0 || phia >= phi(a_lo)
        % Our alpha is an improvement of a_hi
        disp('Case #1');
        a_hi = alpha;
    else
        disp('Case #2');
        phiBara = phiBar(alpha);
        if abs(phiBara) <= -c2*phiBar0   
            disp('Found optimal alpha for c2');
            % success
            break
        end
        if phiBara*(a_hi - a_lo) >= 0
            % Change the interval
            a_hi = a_lo;
        end
        a_lo = alpha;
    end
    c = c + 1;

    % Look out for errors
    if alpha <= 0
        alpha
        c
        phiBar0
        phi0
        a_lo
        a_hi
        error('Zoom: Alpha less than 0')
    end
end

end
