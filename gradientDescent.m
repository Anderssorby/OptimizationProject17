function theta = gradientDescent(l,p)

theta = zeros(1,length(l));
gf = gradf(l,theta,p);
count = 0;
c1 = 0.25;
gf = gradf(l,theta,p);
rho = 0.5;
while count < 1000 && norm(gf)>1E-7
    count = count + 1;
    a=1;
    while f(l,theta-a*gf,p)>f(l,theta,p)-c1*a*(gf*gf')
        a = rho*a;
    end
    theta = theta - gf*a;
    gf=gradf(l,theta,p);
end

count

end
    
