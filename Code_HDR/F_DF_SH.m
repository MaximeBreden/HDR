function [F, DF] = F_DF_SH(u, rho)

% The map F and DF corresponding to steady states of the Swift-Hohenberg
% equation.

N = length(u) - 1; % u = (u_0,...,u_N)

Mu = convo_cos_mat(u); % Multiplication by u operator 
uu = Mu * u; % u^2
uuu = Mu * uu; % u^3
Lap = diag( -(0:N).^2 ); % Laplacian operator
I = eye(N+1); % Identity

F = (-(I+Lap)^2 + rho*I)*u - uuu; %F(u)

Muu = convo_cos_mat(uu); % Multiplication by u^2 operator
DF = (-(I+Lap)^2 + rho*I) - 3*Muu; %DF(u)
