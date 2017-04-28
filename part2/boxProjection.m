function p = boxProjection(g, l, u)
% Used in BCL
if g <= l
    p = l;
elseif g >= u
    p = u;
else
    p = g;
end
end
