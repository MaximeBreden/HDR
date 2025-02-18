function [rmin, rmax] = iscontraction(Y, Z1, Z2, rstar)

%Checking whether the obtained bounds satifsy the assumptions of the
%Newton-Kantorovich theorem which guarantee that we have a contraction

if nargin < 4
    rstar = Inf;
end

rmin = NaN;
rmax = NaN;
if i2f(Z1, 'sup') < 1
    Delta = (1-Z1)^2 - 2*Y*Z2;
    if Delta > 0
        rmin = i2f( (1-Z1-sqrt(Delta)) / Z2 , 'sup');
        rmax = i2f( min( (1-Z1) / Z2, rstar ), 'inf');
        if rmin < rmax
            fprintf("\nSUCCESS :) \n[r_min,r_max) = [%g,%g)\n",rmin,rmax)
        else
            fprintf('\nFAILURE :( \nIt looks like r* is too small\n')
        end
    else 
        fprintf('\nFAILURE :( \n2*Y*Z2 = %f is not smaller than (1-Z1)^2 = %g\n',i2f(2*Y*Z2,'sup'),i2f((1-Z1)^2,'inf'))
    end
else
    fprintf('\nFAILURE :( \nZ1 is larger than 1\n')
end