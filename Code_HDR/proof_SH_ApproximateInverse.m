function [rmin, rmax] = proof_SH_ApproximateInverse(u, rho, nu ,rstar)

% Proof of Theorem 1.3.7, using the Newton-Kantorovich Theorem 1.2.2.

%% Initialization
u = intval(u);
rho = intval(rho);
nu = intval(nu);
rstar = intval(rstar);

N = length(u) - 1;
N3 = 3*N;
N5 = 5*N;

% Padding with zeros to later remove truncation errors
u_N3 = [u;zeros(2*N,1)]; 
u_N5 = [u;zeros(4*N,1)];

% weights for the computation of the \ell^1_\nu norms
weights = [1 2*ones(1,N)] .* nu.^(0:N);
weights_N5 = [1 2*ones(1,N5)] .* nu.^(0:N5);
weights_N3 = weights_N5(1:N3+1);


%% Building A
[~, DF] = F_DF_SH(mid(u), mid(rho)); % This does not need to be rigorous
AN = intval(inv(DF)); % finite part of A
A_tail = 1 ./ (-(1-intval(N+1:N5).^2).^2)'; % tail of A

%% The bound Y
[F_N3, ~] = F_DF_SH(u_N3, rho); % F(\bar{u}). This is rigorous, because the inputs are intervals, and we added enough zeros in u_N3 to prevent truncation errros
AF_N3 = [AN * F_N3(1:N+1); 
         A_tail(1:2*N) .* F_N3(N+2:end)]; % A*F(\bar{u})
Y = weights_N3 * abs(AF_N3)


%% The bound Z1
[~, DF_N5] = F_DF_SH(u_N5, rho); % DF(\bar{u})
ADF_N5N3 = [AN * DF_N5(1:N+1,1:N3+1); 
            diag(A_tail) * DF_N5(N+2:end,1:N3+1)]; % ADF(\bar{u})
B = eye(N5+1, N3+1) - ADF_N5N3; % I - A*DF(\bar{u}). All the columns contributing to the finite part of Z_1
Z1_finite = max( (weights_N5 * abs(B)) ./ weights_N3 );

uu_N3 = convo_cos(u_N3,u_N3);
e_N3 = zeros(N3+1,1);
e_N3(1) = 1;
Z1_tail = weights_N3 * abs(rho*e_N3-3*uu_N3) * abs(A_tail(1));
            
Z1 = max(Z1_finite,Z1_tail)
            
%% The bound Z2
norm_AN = max( (weights * abs(AN) ./ weights) );
norm_A = max(norm_AN, abs(A_tail(1))); % the norm of A
Z2 = 3 * norm_A * ( 2 * weights * abs(u) + rstar )

%% Checking the constraction
[rmin, rmax] = iscontraction(Y, Z1, Z2, rstar);
