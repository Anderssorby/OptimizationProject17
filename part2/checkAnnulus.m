function bool = checkAnnulus(l,p)
%returns true is all points p are inside the annulus spanned by the arms l
bool=1;

lmax = sum(l);
lmin = 2*max(l)-lmax;
for i = 1:length(p(1,:))
    plength = sqrt(p(1,i)^2+p(2,i)^2);
    if plength>lmax || plength<lmin
        bool=0;
        return
    end
end
end