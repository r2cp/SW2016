function [app_lg] = app_log_gamma(z)
%Compute Approximation to log(GAMMA(z))using Stirling's formula
% See Abramowitz and Stegun equation 6.1.37
%Note z can be a vector
%
 const = 1;
 const = const + 1./(12*z);
 const = const + 1./(288*(z.^2));
 const = const - 139./(51840*(z.^3));
 const = const - 571./(2488320*(z.^4));
 app_lg = log(const);
 app_lg = app_lg+log(sqrt(2*pi));
 app_lg = app_lg + (z-0.5).*log(z);
 app_lg = app_lg -z;
end

