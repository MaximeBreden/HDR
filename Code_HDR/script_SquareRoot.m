clear variables
close all
clc
format long

%% Computation of an approximation x of sqrt(2)
F = @(x) x.^2 - 2;
DF = @(x) 2*x;

x = 1;
err = abs(F(x));
it = 0;
tol = 1e-15;
while err > tol && it < 10 
    x = x - DF(x)\F(x);
    err = abs(F(x));
    it = it+1;
end

x

%% Validation of the approximation using the Newton-Kantorovich theorem
x = intval(x);
A = 1/(2*x);
Y = abs(A*F(x));
Z1 = 0;
Z2 = 1/abs(x);

fprintf('\n Y = %.7e\n', sup(Y))
fprintf('\n Z2 = %.7e\n', sup(Z2))
[rmin, rmax] = iscontraction(Y, Z1, Z2);

