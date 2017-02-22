function alpha = zoom(phi, phiBar, constants, a_lo, a_hi)
c1 = constants(1);
c2 = constants(2);

phiBar0 = phiBar(0);
phi0 = phi(0);

c = 0;

while c < 1000
    % Interpolate to find a test step length alpha

    % Minimizing the interpolation on [a_lo, a_hi]
    alpha = (phiBar(a_lo)*a_hi^2)/(2*(phi(a_hi) - phi(a_lo)- phiBar(a_lo)*a_hi));


    phia = phi(alpha);
    
    if phia > phi0 + c1*alpha*phiBar0 || phia >= phi(a_hi)
        a_hi = alpha;
    else
        phiBara = phiBar(alpha);
        if abs(phiBara) <= -c2*phiBar0   
            % success
            break
        end
        if phiBara*(a_hi-a_lo) >= 0
            a_hi = a_lo;
        end
        a_lo = alpha;
    end
    c = c + 1;
end

end
