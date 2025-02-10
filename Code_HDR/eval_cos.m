function ux = eval_cos(u, x)

N = length(u) - 1;

if size(x,1) == 1
    x = transpose(x);
end

u(2:end) = 2 * u(2:end);
ux = cos(x*(0:N)) * u;
