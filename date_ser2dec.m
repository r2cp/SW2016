function dec_date = date_ser2dec(serdate);

% This function converts a Matlab date, in serial number form, to a decimal date, for example 1971.246.
% Some notes
% ------------ input ------------
% serdate = Matlab serial date vector, Jan 1, 0000 AD = 1
%           Note Excel (Mac 2011) and Windows use Jan 1, 1900 = 1
%                Excel (Mac 2008) uses Jan 1, 1904 = 1
% ------------ output ------------
% dec_date = decimal date (e.g., 1971.246)
%
% Find serial number for first year and last year
dvec = datevec(min(serdate));
fyear = dvec(1,1);
dvec = datevec(max(serdate));
lyear = dvec(1,1);
first_ser = datenum(fyear,1,1);
% Use this as new base
serdate = serdate - first_ser + 1;
reg_year = [1:1:365]'/365;
leap_year = [1:1:366]'/366;
dec_vec = 0;
for iyear = fyear:1:lyear;
	i_leap = leap_check(iyear);
  if i_leap == 0;
   tmp = reg_year;
  elseif i_leap == 1;
   tmp = leap_year;
  end;
  tmp = tmp+iyear;
  dec_vec = [dec_vec;tmp];
end;
dec_vec = dec_vec(2:end);
dec_date = dec_vec(serdate);