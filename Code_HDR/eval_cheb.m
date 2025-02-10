function xlambda = eval_cheb(x, lambda)

K = length(x) - 1;

if size(lambda,1) == 1
    lambda = transpose(lambda);
end

x(2:end) = 2 * x(2:end);
xlambda = cos(acos(lambda)*(0:K)) * x;
