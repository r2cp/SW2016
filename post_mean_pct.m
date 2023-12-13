function [x_output] = post_mean_pct(x,pctvec)
% Compute posterior mean and percentiles for the columns of x
  %  x is (n_draws) x (n_s), where n_draws is the number of draws
  %  and n_s is the number of series
  %  x_output is the output vector
  %   col(1) = mean
  %   cols(2)-end = percentiles
  %   number of rows = n_s
  n_s = size(x,2);
  n_p = size(pctvec,1);
  x_output = NaN*ones(n_s,n_p+1);
  for i = 1:n_s;
      tmp = x(:,i);
      x_output(i,1) = mean(tmp);
      pct_tmp = pctile(tmp,pctvec);
      x_output(i,2:end) = pct_tmp';
  end;

end

