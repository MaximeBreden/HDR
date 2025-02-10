function kappa0 = compute_kappa0(u, rho, show_enclosures)

% Computation of kappa0, using the homotopy method

N = length(u) - 1;

op = initialize_L(u, rho); % Construct the operators L and L0, shifted so that all the eigenvalues are positive, and a couple of auxiliairy things)

%% Determing how many eigenvalues of the base problem to start with
M0 = 1;
M0largeenough = false;
while not(M0largeenough) && M0 <= N+1
    eigenval = eigs(op.L, M0, 'SM'); % first eigenvalues of L + shift
    eigenval = sort(real(eigenval)) - op.shift; % first eigenvalues of L
    M0largeenough = (eigenval(1) > 0) || (eigenval(M0) > 0) ; % We want to have enough eigenvalues of L to have those around zero
    if not(M0largeenough)
        M0 = M0+1;
    end
end

if not(M0largeenough)
    fprintf("No suitable M found")
    kappa0 = Inf;
    return
end

target_ev = eigenval(M0); % The first eigenvalue of L that is larger than 0.
margin = 1; % We are going to take enough eigenvalues of L0 so that the largest is at greater than target_ev + margin
mmax = ceil(sqrt(1+sqrt(1+op.cmax+target_ev+margin))); % We should only need the eigenvalues of L0 for 0 <= m <= mmax
m = intval(0:mmax);
ieigenval_0 = m.^4 - 2*m.^2 - op.icmax + op.ishift; % rigorous eigenvalues of L0 + shift
eigenvect_0 = eye(N+1,mmax+1); % corresponding eigenvectors
for m = 2 : mmax+1 % normalization of the eigenvectors
    eigenvect_0(m,m) = 1/sqrt(2);
end

% Sorting the eigenvalues (and associate eigenvectors), so that the
% eigenvalue are in ascending order
[~, ind] = sort(mid(ieigenval_0));
ieigenval_0 = ieigenval_0(ind);
eigenval_0 = mid(ieigenval_0);
eigenvect_0 = eigenvect_0(:,ind);


%% Using the homotopy method to rigorously enclose eigenvalues of L
s = 0; % starting point for the homotopy
steps_s = 0;
margin_L = 0.01; %This parameter controls how the breakpoints s_k for the homotopy are chosen.
%We take s_{k+1} such that \lambda^{(s_{k+1})}_{M-(k+1)} \approx (1-margin_L) * \lambda^{(s_{k})}_{M-k}.
%This may have to be chosen more carefully for more complicated problems.

para_homotopy.tol_simple_eig = 10^-7; %For the validation of simple VS clustered eigenvalues
para_homotopy.show_s = 0; %Put to 1 to display what the value of s is during the whole process
para_homotopy.show_enclosure = show_enclosures; %Put to 1 to display the obtained enclosures at each homotopy breakpoint

% The whole homotopy method (including the Rayleigh-Ritz and Lehmann-Maehly
% methods), happens in the function "homotopy" called below.
[ilower_eig_L,iupper_eig_L,~,~] = homotopy(s, steps_s, op, eigenval_0, eigenvect_0, ieigenval_0.inf, margin_L, para_homotopy);


eig_L = infsup(ilower_eig_L, iupper_eig_L) - op.ishift; %rigorous enclosures of the first eigenvalues of L
fprintf('\nRigorous enclosure of the first eigenvalues of L:\n')
infsup(eig_L)

if inf(eig_L(1)) > 0
    kappa0 = 1/eig_L(1);
elseif max(inf(eig_L)) > 0
    kappa0 = 1/min(abs(eig_L));
else
    fprintf("We don not have enough eigenvalues of L to control kappa0")
    kappa0 = Inf;
end
    
    

