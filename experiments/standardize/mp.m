function [ density ] = mp( x, gamma )

% Analytical Marcenko-Pastur density
a = @(gamma) ((1-sqrt(gamma))^2);
b = @(gamma) ((1+sqrt(gamma))^2);
density= (1./(2*pi*x*gamma)).*sqrt((b(gamma)-x).*(x-a(gamma)));

end

