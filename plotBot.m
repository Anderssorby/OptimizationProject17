function plotBot(l,theta)

cords = zeros(2,length(l)+1);


figure; hold on;

for i = 1:length(l)
    cords(:,i+1)=transpose(bigF(l,theta,i));
end

for i = 1:length(l)
    plot([cords(1,i), cords(1,i+1)],[cords(2,i), cords(2,i+1)],'bl')
    plot(cords(1,i),cords(2,i),'k*')
end
plot(cords(1,end),cords(2,end),'k*')
grid on
xlabel('x'), ylabel('y')
title('Robot arm')
end
