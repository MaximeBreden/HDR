function xcheb = interp_cheb(xnodes)

K = length(xnodes) - 1;

if isa(xnodes(1), 'interval')
    ipi = intval('pi');
else
    ipi = pi;
end

theta_K = (K-(0:K)')*ipi/K;

% Construction of the matrix to go from values at Chebyshev nodes to
% Chebyshev coefficients.
% For this simple problem, using the FFT is not needed, and instead we
% simply build the matrix corresponding to the DCT, which allows to pass
% from values at Chebyshev nodes to Chebyshev coefficients.
Minterp = transpose(cos(theta_K*(0:K))) / K; 
Minterp(:,[1,K+1]) = Minterp(:,[1,K+1]) / 2;
Minterp(K+1,:) = Minterp(K+1,:) / 2;

xcheb = Minterp * xnodes;