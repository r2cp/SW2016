function xbp = bpassar_1sd(x,updc,lpdc,n,nar,tcode,nobs_min)
%   1-sided bandpass ... constructed using repeated calls to bpassar
%   Input:
%     x     = series to be filtered
%     updc  = period corresponding to upper cutoff frequency
%     lpdc  = period corresponding to lower cutoff frequency
%     n     = numbers of terms in moving average filter
%     nar   = order of AR used for padding
%     tcode = transformation code for AR model
%              0 -- no transformation
%              1 -- first difference
%     nobs_min = minimum number of obs for compuatiion
%
%   Output:
%    xfiltered = filtered value of series, x
dnobs = size(x,1);
xbp = NaN*zeros(dnobs,1);
for t = nobs_min:dnobs;
 xp = packr(x(1:t));
 if size(xp,1) >= nobs_min;
   [xbp_t, sebp_t ] = bpassar(x(1:t),updc,lpdc,n,nar,tcode,0);
   xbp(t)=xbp_t(end);
 end;
end