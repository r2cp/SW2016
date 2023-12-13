function [covmat] = acv2cov(acv)
% Construct a stationary covariance matrix from a vector of autovariances
%
 r = size(acv,1);
 covmat = NaN*zeros(r,r);
 covmat(:,1)=acv;
 covmat(:,end)=acv(end:-1:1);
 for i = 2:r-1;
    covmat(:,i)=[acv(i:-1:1);acv(2:(r-i+1))];
 end;

end

