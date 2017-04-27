function a1 = interpolateQuad(a, b, f_a, fB_a, f_b)
% Interpolating the points f(a), f'(a) and f(b) 

a1 = (fB_a*b^2)/(2*(f_b - f_a - fB_a*b));

end
