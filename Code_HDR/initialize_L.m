function Lstruct = initialize_L(u, rho)

%% Construction of the basic objects involved in the definition of L and L0
N = length(u) - 1;

iu = [intval(u); zeros(2*N,1)]; % Adding enough zero to remove all truncations errors
irho = intval(rho);

uu = convo_cos(u,u); % u^2
e = zeros(N+1,1);
e(1) = 1;
c = (rho - 1)*e - 3*uu;
Mc = convo_cos_mat(c); % Multiplication by c operator
Lap = diag( -(0:N).^2 ); % Laplacian operator
I = eye(N+1); % Identity

iuu = convo_cos(iu,iu); % u^2 
ie = zeros(3*N+1,1);
ie(1) = 1;
ic = (irho - 1)*ie - 3*iuu;
iMc = convo_cos_mat(ic); % Multiplication by c operator
iLap = diag( -(0:3*N).^2 ); % Laplacian operator
iI = eye(3*N+1); % Identity
ipi = intval('pi');

% Rigorous bound cmax >= max_{x in [0,pi]} u(x). 
x = linspace(0, pi, 1e4);
x(1) = 0;
x(end) = sup(ipi);
ix = infsup(x(1:end-1), x(2:end)); % Subdivision of [0,pi] into smaller subintervals
cmax = sup(max(eval_cos(ic, ix))); % This is not very efficient (but more than good enough for this simple problem) 
icmax = intval(cmax);


%% We shift all the spectra to be safely away from zero. 
% The parameter shift is chosen so that all the eigenvalues of L^0 + shift are >= 1.
shift = cmax + 2;
ishift = icmax + 2;

%% Construction of L and L0 (with extented version to remove truncation errors)
L = Lap^2 + 2*Lap - Mc + shift*I;
L0 = Lap^2 + 2*Lap - cmax*I + shift*I;

% Whenever v is a vector of size N (padded with 2N extra zeros), Lv and Lv0
% are computed exactly by doing iL*v and iL0*v. 
iL = iLap^2 + 2*iLap - iMc + ishift*iI;
iL0 = iLap^2 + 2*iLap - icmax*iI + ishift*iI;


% Matrices of the scalar product 
Mat_scal = diag([1 2*ones(1,N)]);
iMat_scal = diag([1 2*ones(1,3*N)]);


%% Structure containing everything that will be needed for the homotopy method
Lstruct.L = L;
Lstruct.L0 = L0;
Lstruct.Mat_scal = Mat_scal;
Lstruct.Id = I;
Lstruct.shift = shift;
Lstruct.cmax = cmax;

Lstruct.iL = iL;
Lstruct.iL0 = iL0;
Lstruct.iMat_scal = iMat_scal;
Lstruct.ishift = ishift;
Lstruct.icmax = icmax;
