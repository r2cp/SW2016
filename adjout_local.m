function [x] = adjout_local(y,thr,tflag,bwidth)
% -- Adjust for outliers using fraction of IQR
%      
%       y = Data series
%       thr = threshold in multiples of IQR
%       tflag = 0  == replace with missing value 
%               1  == replace with median value in rolling window
%       bwidth: data are checked over rolling periods of this bandwidth
%

small = 1.0e-06;

% -- Compute IQR 
pct_vec = [0.25 0.50 0.75];
bw2 = floor(bwidth/2);
x = NaN(size(y,1),1);

for t = 1:size(y);
    if t < 1+bw2;
        zt = y(1:bwidth);
    elseif t > size(y,1)-bw2;
        zt = y(size(y,1)-bwidth+1:end);
    else;
        zt = y(t-bw2:t+bw2);
    end;
    z = zt(isnan(zt) == 0);
    tmp = pctile(z,pct_vec);
    zm = tmp(2);
    iqr = tmp(3)-tmp(1);     
    if iqr < small;
      x(t) = NaN;
    end;
    if iqr >= small;
      x(t) = y(t);
      ya=abs(y(t)-zm);
      iya = ya > (thr*iqr);
      if iya == 1;
        if tflag == 0;
           x(t) = NaN;
        elseif tflag == 1;
           x(t)=zm;
        end;
      end;
    end;

end;


end

