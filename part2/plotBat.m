function plotBat(l,theta,pmat)

cords = zeros(2,length(l)+1);


figure; hold on;
for j = 1:length(pmat(1,:))
    
    for i = 1:length(l)
        cords(:,i+1)=transpose(bigEff(l,theta((j-1)*length(l)+1:j*length(l)),i));
    end

    for i = 1:length(l)
        plot([cords(1,i), cords(1,i+1)],[cords(2,i), cords(2,i+1)],'bl')
        plot(cords(1,i),cords(2,i),'ko')
    end
    plot(cords(1,end),cords(2,end),'k')
    %grid on
end
    
    

xlabel('x'), ylabel('y')
title('Robot arm')
for i = 1:length(pmat(1,:))
    plot(pmat(1,i),pmat(2,i),'r*')
end
axis equal
end
    
