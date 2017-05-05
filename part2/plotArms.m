function plotArms(l,theta,p,holdon)
%Plots all arms and points for a given l, theta and p. 
%holdon is a bool telling if this should be plotted in a figure that is
%already open.

if ~holdon
    figure; 
end
hold on;   

cords = zeros(2,length(l)+1); 

for j = 1:length(p(1,:))
    
    for i = 1:length(l)
        cords(:,i+1)=transpose(bigEff(l,theta((j-1)*length(l)+1:j*length(l)),i));
    end

    for i = 1:length(l)
        plot([cords(1,i), cords(1,i+1)],[cords(2,i), cords(2,i+1)],'black','LineWidth',2)
        plot(cords(1,i),cords(2,i),'ko')
    end
    plot(cords(1,end),cords(2,end),'k')
    %grid on
end
    


xlabel('x'), ylabel('y')
title('Robot arm')
for i = 1:length(p(1,:))
    plot(p(1,i),p(2,i),'r*')
end
axis equal
end
    
