
% column vectors
fprintf('#1 Test \n')
l = [3; 2; 2]
p = [3;2]

x0 = [1; 1; 2];
[theta, k] = BFGS(l, p, x0);
thGD = gradientDescent(l, p, x0);
fprintf('Solution BFGS\n');
theta
f(l,theta',p)
fprintf('Solution gradient descent\n');
thGD
f(l,thGD, p)



fprintf('#2 Testing l = %s, p = %s\n',l,p)
l = [1; 4; 1]
p = [1;1]

x0 = [1; 1; 2];
[theta, k] = BFGS(l, p, x0);
thGD = gradientDescent(l, p, x0);
fprintf('Solution BFGS\n');
theta
f(l,theta',p)
fprintf('Solution gradient descent\n');
thGD
f(l,thGD, p)
%k
%p
%bigF(l,theta)


l = [3; 2; 1; 1];
p = [3;2];

x0 = [1; 1; 2; 0];
fprintf('#3 Testing\n')
[theta, k] = BFGS(l, p, x0);
thGD = gradientDescent(l, p, x0);
fprintf('Solution BFGS\n');
theta
f(l,theta',p)
fprintf('Solution gradient descent\n');
thGD
f(l,thGD, p)
%k
%p
%bigF(l,theta)


fprintf('#4 Test')
l = [3; 2; 1; 1]
p = [0;0]

x0 = [1; 1; 2; 0];
[theta, k] = BFGS(l, p, x0);
thGD = gradientDescent(l, p, x0);
fprintf('Solution BFGS\n');
theta
f(l,theta',p)
fprintf('Solution gradient descent\n');
thGD
f(l,thGD, p)
%k
%p
%bigF(l,theta)
