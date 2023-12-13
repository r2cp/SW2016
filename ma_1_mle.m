function [theta_hat,sigma2_e] = ma_1_mle(y,theta_min,theta_max,grid_value);
% Compute MLE for MA(1) Model
% y(t) = e(t) - theta*e(t-1)
% MLE constructed by grid search over values of theta
% theta_min = minimum value of theta
% theta_max = maximum value of theta
% grid_value = length of grid for theta search

theta_grid = (theta_min:grid_value:theta_max);
n_grid = size(theta_grid,2);
T = size(y,1);
e_std = NaN*ones(T,n_grid);
for i = 1:n_grid;
    theta = theta_grid(i);
    H = [1 -theta]';
    F = zeros(2,2);
    F(2,1) = 1;
    Q = zeros(2,2);
    Q(1,1) = 1;
    p1 = eye(2);
    p2 = zeros(2,2);
    p2(1,1) = 1;
    x1 = zeros(2,1);
    x2 = zeros(2,1);
    for t = 1:T;
        x2 = F*x1;
        p2 = F*p1*F' + Q;
        h = H'*p2*H;
        e = y(t)-H'*x2;
        k = p2*H/h;
        x1 = x2 + k*e;
        p1 = (eye(2) - k*H')*p2;
        e_std(t,i) = e/sqrt(h);
     end;
end;
ssr = sum(e_std.^2)';
[ssr_min,ii] = min(ssr);
theta_hat = theta_grid(ii);
sigma2_e = ssr_min/T;

end

       
        



