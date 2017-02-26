% a = 100;
% l = randi(20,a,1);
% p = [randi(20), randi(30)]';
% x0 = 2*pi*rand(a,1);

l = [1 3 2 1]';
p = [0 1]';
x0 = [1 1 1 1]';

[u,vbar,wbar] = gradientDescent(l,p,x0);
[U, Vbar, Wbar]=BFGS2(l,p,x0);
for i = 1:10
[u,v,w] = gradientDescent(l,p,x0);
[U, V, W]=BFGS2(l,p,x0);
vbar = (i*vbar+v)/(i+1);
Vbar = (i*Vbar+V)/(i+1);
end

n = length(v);
N = length(V);

figure;
subplot(1,2,1);
semilogx(w,1:n);
hold on
semilogx(W,1:N);
title('Error vs. steps')
axis tight
legend('Gradient descent','BFGS')
ylabel('Number of steps')
xlabel('Error')

subplot(1,2,2);
semilogx(w,vbar);
hold on
semilogx(W,Vbar);
ylabel('Time spent')
xlabel('Error')
title('Error vs. time')
axis tight
legend('Gradient descent','BFGS')




