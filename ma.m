function [xma] = ma(x,p)
% Form p term MA of X
 xma = x;
 for i=1:p-1;
 	xma = xma+lag(x,i);
 end;
 xma = xma/p;
end