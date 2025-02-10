clear variables
close all
clc
format long

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% You need to use Intlab to run this script %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% The map f
f = @(x) exp(sin(20*x)) + x.^4 + x./(1+x.^2) - 1/8;

%% Initial plot (nonrigorous)
x_plot = linspace(-1, 1, 1e3);
figure
plot(x_plot, 0+0*x_plot, 'k')
hold on
plot(x_plot, f(x_plot), 'b', 'linewidth', 2)
xlabel('$x$', 'Interpreter', 'Latex')
title('Graph of $f$', 'Interpreter', 'Latex')
set(gca, 'FontSize', 15) 
axis tight

%% First rigorous enclosures
fprintf("\n<strong>Example 1.1.1.</strong>")
x0 = intval('-0.4');% Creates an interval containing the real number -0.4,
% which is not a representable number because the computer uses base 2.
fprintf("\nGuaranteed enclosure of f(-0.4):\n")
infsup(f(x0))


fprintf("\n<strong>Remark 1.1.2.</strong>")
I = infsup(0.25,0.5);
fprintf("\nGuaranteed enclosure of f([0.25,0.5]):\n")
f(I)

fprintf("\n<strong>Example 1.1.3.</strong>")
x1 = intval('-0.45');
fprintf("\nGuaranteed enclosure of f(-0.45):\n")
infsup(f(x1))

x2 = intval('-0.35');
fprintf("\nGuaranteed enclosure of f(-0.35):\n")
infsup(f(x2))


%% Sharper enclosures
fprintf("\n<strong>Remark 1.1.4.</strong>")
% One could do many different things here, e.g., use Interval Newton). We
% simply compute approximate zeros using the standard Newton method, and
% then slightly perturb them to the left and to the right.

df = @(x) 20*exp(sin(20*x)).*cos(20*x) + 4*x.^3 + (1-x.^2)./(1+x.^2).^2; % The derivative of f
x1float = -0.45;
x2float = -0.35;
for k = 1:10 % A few interations of Newton's method on both approximate zeros
    x1float = x1float  - f(x1float)/df(x1float);
    x2float = x2float  - f(x2float)/df(x2float);
end

x1left = intval(x1float - eps); % Should be to the left of the exact zero
x1right = intval(x1float + eps); % Should be to the right of the exact zero
if f(x1left)*f(x1right) < 0 % We check that the intermediate value theorem applies
    x1sharp = infsup(inf(x1left),sup(x1right));
    fprintf("\nSharper interval containing a zero of f:\n")
    infsup(x1sharp)
end

x2left = intval(x2float - eps); % Should be to the left of the exact zero
x2right = intval(x2float + eps); % Should be to the right of the exact zero
if f(x2left)*f(x2right) < 0 % We check that the intermediate value theorem applies
    x2sharp = infsup(inf(x2left),sup(x2right));
    fprintf("\nSharper interval containing another zero of f:\n")
    infsup(x2sharp)
end


%% Global result
fprintf("\n<strong>Example 1.1.5.</strong>")
N = 2^8; % Number of subintervals in the subdivision
grid = linspace(-1, 1, N+1);
intervals = infsup(grid(1:end-1), grid(2:end)); % The subintervals
fintervals = f(intervals); % Enclosures of f on each subintervals

figure
plot([-1 1], [0 0], 'k')
hold on
plot(intervals, fintervals)
xlabel('$x$', 'Interpreter', 'Latex')
title('Guaranteed enclosures of $f$', 'Interpreter', 'Latex')
set(gca, 'FontSize', 15) 
axis tight

sign_definite = max(fintervals > 0, fintervals < 0); % Checking which intervals are guaranteed to not contain a zero of f
potential_zeros = find(sign_definite == false); % Intervals which may still contain zero(s) of f
nb_int = length(potential_zeros);
fprintf("\nThere are at most %d intervals containing zeros of f:\n", nb_int)
for k = 1:nb_int
    fprintf("\nThe %d-th one:\n", potential_zeros(k))
    infsup(intervals(potential_zeros(k)))
end

% min( f(intervals(1:73)) > 0 )
% min( f(intervals(75:81)) < 0 )
% min( f(intervals(83:end)) > 0 )

dfintervals = df(intervals(potential_zeros));
if max(dfintervals > 0, dfintervals < 0) % Checking that f is monotone on each remaining subinterval
    fprintf("\nf is monotone on each of these %d intervals.\n", nb_int)
end



