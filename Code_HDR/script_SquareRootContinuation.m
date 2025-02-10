clear variables
close all
clc
format short

%% Approximate zero
% We look for a function x(lambda) so that x(lambda)^2 = 2 + lambda,
% for all lambda in [-1,1].

K = 20; % number of Chebyshev modes used for the approximate solution
nodesK = cos((K:-1:0)*pi/K);
x_atthenodes = zeros(K+1,1);

DF = @(x) 2*x;
tol = 1e-15;
for k = 1 : K+1
    lambda = nodesK(k);
    F = @(x) x.^2 - (2 + lambda);
    x = 1;
    err = abs(F(x));
    it = 0;
    while err > tol && it < 10 
        x = x - DF(x)\F(x);
        err = abs(F(x));
        it = it+1;
    end
    x_atthenodes(k) = x;
end

% the approximate zero of F in ell^1_\eta
xbar = interp_cheb(x_atthenodes);

% Sanity check
lambda = linspace(-1, 1, 1e3);
figure
plot(2+lambda, eval_cheb(xbar,lambda) ,'b', 'Linewidth', 2)
hold on
plot(2+lambda, sqrt(2+lambda) ,'r--', 'Linewidth', 2)
xlabel('$\lambda$', 'Interpreter', 'Latex')
legend('$\bar{x}(\lambda)$', '$\sqrt{\lambda}$', 'Interpreter', 'Latex', 'Location', 'SouthEast')
title('Approximation of the square root function')
set(gca, 'FontSize', 15) 
axis tight

%% Proof
% the ybar such that A = M(ybar)
ybar = interp_cheb(1./(2*x_atthenodes));

% Padding by zeros to remove truncation errors
xbar = [intval(xbar); zeros(2*K,1)];
ybar = [intval(ybar); zeros(2*K,1)];
g = zeros(3*K+1,1);
g(1) = 2;
g(2) = 1/2;

% weights for the norms
eta = intval('1'); % weight for the \ell^1_\eta space
weights = [1 2*eta.^(1:3*K)];

% Y bound
nodes3K = cos((3*K:-1:0)*intval('pi')/(3*K));
AFx_atthenodes = eval_cheb(ybar,nodes3K) .* (eval_cheb(xbar,nodes3K).^2 - eval_cheb(g,nodes3K));
AFx = interp_cheb(AFx_atthenodes);
Y = weights * abs(AFx)

% Z1 bound (b = 1-2x*y)
e = zeros(3*K+1,1);
e(1) = 1;
b_atthenodes = eval_cheb(e,nodes3K) - 2 * eval_cheb(ybar,nodes3K) .* eval_cheb(xbar,nodes3K);
b = interp_cheb(b_atthenodes);
Z1 = weights * abs(b)

% Z2 bound
Z2 = 2 * weights * abs(ybar)

% Checking the assumptions of the Newton-Kantorovich theorem
[rmin, rmax] = iscontraction(Y, Z1, Z2);








