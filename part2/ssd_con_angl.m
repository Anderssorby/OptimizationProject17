function E = ssd_con_angl(Theta)
% Sums of squares of distances between consecutive angles
[n, s] = size(Theta);

for i=1:n
   for j=1:s-1
      E = E + (Theta(i, j+1) - Theta(i, j))^2;
   end
end

for i=1:n
    E = E + (Theta(i, 1) - Theta(i, s))^2;
end

E = 0.5*E;
end
