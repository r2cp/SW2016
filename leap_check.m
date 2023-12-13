function i_leap = leap_check(iyear);

% Check to see if iyear is a leap year

small = .00000001;
i_4 = 0;
i_100 = 0;
i_400 = 0;
tmp = iyear/4;
tmp1 = tmp-floor(tmp);
if tmp1 < small;
  i_4 = 1;
end;
tmp = iyear/100;
tmp1 = tmp-floor(tmp);
if tmp1 < small;
  i_100 = 1;
end;
tmp = iyear/400;
tmp1 = tmp-floor(tmp);
if tmp1 < small;
  i_400 = 1;
end;
i_leap = i_4*(1-i_100)+i_400;
