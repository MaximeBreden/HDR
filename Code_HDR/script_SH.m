clear variables
close all
clc
format short

%% INITIALIZATION
load('data_SH.mat','u','rho') % A precomputed approximate solution

% Parameters for the proof
nu = intval('1');
rstar = intval('1e-5');

% Plot
x = linspace(0, pi, 1e3)';
figure
plot(x, eval_cos(u,x), 'Linewidth', 2)
legend('$\bar{u}(x)$', 'Interpreter', 'Latex', 'Location', 'NorthWest')
xlabel('$x$', 'Interpreter', 'Latex')
title(['Approximate solution for $\rho$ = ',num2str(rho)], 'Interpreter', 'Latex')
set(gca, 'FontSize', 15) 
axis tight
drawnow

%% PROOF WITH APPROXIMATE INVERSE
fprintf("\nFirst proof with an approximate inverse, computation of the bounds ...\n")
[rmin1, rmax1] = proof_SH_ApproximateInverse(u, rho, nu, rstar);

%% PROOF WITH EIGENVALUE BOUNDS
fprintf("\n\nSecond proof with eigenvalue bounds, computation of the bounds ...\n")
show = 0; % Change to show = 1 in order to get intermediate information during the homotpy method
[rmin2, rmax2] = proof_SH_EigenvalueBounds(u, rho, rstar, show);
