function confSpace(l,MA,holdon)
n = length(l);
nAngles = 30;
angles = -MA:2*MA/(nAngles-1):MA;

theta = zeros(n,1);
pmat = zeros(2,nAngles^n);

countvec = ones(1,n);

j=0;
while sum(countvec) < n*nAngles
    j = j+1;
    for i = 1:n
        theta(i) = angles(countvec(i));
    end
    
    pmat(:,j) = bigEff(l,theta,n);
    
    k = n;
    while countvec(k) == nAngles
        countvec(k) = 1;
        k = k-1;
    end
    countvec(k) = countvec(k) + 1;
    
end
if~holdon
    figure;
end
plot(pmat(1,:),pmat(2,:),'.');
axis equal
grid on
hold on

end