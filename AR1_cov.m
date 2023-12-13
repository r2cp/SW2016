function [covmat] = AR1_cov(T,rho)
%  Compute TxT covariance matrix for stationary AR(1) with coef = rho and
%  unit variance of innovation
%
 acv = rho.^(0:T-1)';
 acv = acv/(1-rho^2);
 covmat = acv2cov(acv);
end

