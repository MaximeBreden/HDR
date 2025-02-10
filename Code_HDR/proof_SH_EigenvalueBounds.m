function [rmin, rmax] = proof_SH_EigenvalueBounds(u, rho, rstar, show)

% Proof of Theorem 1.3.18, using the Newton-Kantorovich Theorem 1.2.6.

N = length(u) - 1;
iu = intval(u);
ipi = intval('pi');

%% The bound kappa
% Computation of kappa0 using the homotopy method
kappa0 = compute_kappa0(u, rho, show)

% Preliminaries for the computation of kappa1
iu_N2 = [iu; zeros(N,1)];
irho = intval(rho);
iuu = convo_cos(iu_N2,iu_N2); % \bar{u}^2
ie = zeros(2*N+1,1);
ie(1) = 1;
ic = (irho - 1)*ie - 3*iuu; % c = rho - 1 - 3\bar{u}^2

x = linspace(0, pi, 1e4);
x(1) = 0;
x(end) = sup(ipi);
ix = infsup(x(1:end-1), x(2:end)); % Subdivision of [0,pi] into smaller subintervals

% Computation of a rigorous upper bound for || c ||_\infty
norminf_c = max(abs(eval_cos(ic, ix))); % This is not very efficient (but more than good enough for this simple problem) 

% The function defining kappa1 in terms of kappa0
func = @(eta) sqrt( ( 1 + typeadj(kappa0,typeof(eta))^2*( typeadj(norminf_c,typeof(eta))^2*(1/eta(2)+1/(2*eta(3))) + 1/(2*eta(4))*(2/eta(1)+2*eta(3)) ) ) ...
                  / ( 1 - (2*eta(1)+eta(2)+eta(4)/2) ) );
constraint = @(eta) 2*eta(1) + eta(2) + eta(4)/2 * (2/eta(1) + 2*eta(3)) - 1;

% Some suitable coefficients eta
eta = [1/16; 1/8; 16; 1/128];

% Optimizaton to find better values of eta
options = optimoptions(@fmincon,'Display','off');
eta = optimize_eta(eta); 

eta = intval(eta);
if min(inf(eta)) <= 0 || sup(constraint(eta)) >= 0
    fprintf("Invalid choice of eta")
    rmin = NaN;
    rmax = NaN;
    return
end
 
% Computation of kappa1 and kappa
kappa1 = func(eta);
kappa = sqrt(kappa0^2+kappa1^2)

%% The bound delta
iu_N3 = [iu; zeros(2*N,1)];
irho = intval(rho);
[F_N3, ~] = F_DF_SH(iu_N3, irho); % F(\bar{u})
delta = sqrt(F_N3(1).^2 + 2*sum(F_N3(2:end).^2)) % || F(\bar{u}) ||_{L^2}

%% The bound L

normu = sqrt(iu(1).^2 + 2*sum(iu(2:end).^2)); % || \bar{u} ||_{L^2}
L = 6 * (ipi^4 + ipi^2/2)^2 * ( 2*normu + rstar)

%% Checking the constraction
Y = delta * kappa;
Z1 = 0;
Z2 = L * kappa;
[rmin, rmax] = iscontraction(Y, Z1, Z2, rstar);
                                            
              
function  out = optimize_eta(eta0)    
    out = fmincon(func, eta0, [], [], [], [], zeros(4,1), Inf(4,1), @nonlcon, options);
    
    function [cineq, ceq] = nonlcon(eta)
        cineq = constraint(eta);
        ceq = [];
    end
end
    

end