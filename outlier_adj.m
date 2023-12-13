function x_adj = outlier_adj(x,bw,outlier_mult,repl_value);
% Inputs:
%          x = series
%          bw = bandwidth
%          outlier_mult = number of IQRs for outlier detection
%          repl_value = 1 (median over BW period)
%                       0  Missing value
%
 n = size(x,1);
 x_adj = NaN*ones(n,1);
 bw = 41;  % bandwidth;
 bw2 = floor((bw-1)/2);
 for t = 1:n;
 	t1 = max([1,t-bw2]);
 	t2 = min([n,t+bw2]);
 	z = x(t1:t2);
 	zpct = pctile(z,[0.25,0.50,0.75]);
 	zmed = zpct(2);
 	ziqr = zpct(3)-zpct(1);
 	xnorm = (x(t)-zmed)/ziqr;
 	if abs(xnorm) < outlier_mult;
 	  x_adj(t) = x(t);
 	else;
 	  if repl_value == 1;
 	   x_adj(t) = zmed;
 	  else;
 	   x_adj(t) = NaN;
 	  end;
 	end;
 end;

end