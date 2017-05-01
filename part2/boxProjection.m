function g = boxProjection(g, l, u)
% Used in BCL
for i = 1:length(g)
    if g(i) <= l
        g(i) = l;
    elseif g(i) >= u
        g(i) = u;
    else
        g(i) = g(i);
    end
end