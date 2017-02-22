function theta = gradientDescent(l,p)

theta = zeros(1,length(l));
count = 0;
c1 = 0.01;
gf = gradf(l,theta,p);
rho = 0.5;
while 1
    count = count + 1;
    a=1;
    while f(l,theta-a*gf,p)>f(l,theta,p)-c1*a*(gf*gf')
        a = rho*a;
    end
    theta = theta - gf*a;
    gf=gradf(l,theta,p);
    if count > 1000
        error('Surpassed max number of iterations');
    end
    if norm(gf)<1E-7
        [b,theta] = checkEnd(l,theta,p);
        if b
            break
        end
    end
end
end
