function g = boxProjection(g, l, u)
% Projects all the values in a vector onto the interval decided by l(lower)
% and u(upper)
for i = 1:length(g)
    if g(i) <= l
        g(i) = l;
    elseif g(i) >= u
        g(i) = u;
    else
        g(i) = g(i);
    end
end