% Set up data and calendars for monthly data in pcomp .. 
% Construct Monthly and Quarterly Values of Inflation
% 9/9/2015, mww

if load_data == 1;

% ----------- Features of Data Set ---------
% ns = 26;    % Number of PCE series in the files
% miss_code = 1.0e+32;

% ----------- Sample Period, Calendars and so forth
[dnobs_m,calvec_m,calds_m] = calendar_make([2000 12], [2023 9], 12);
[dnobs_q,calvec_q,calds_q] = calendar_make([2000 12] , [2023 9], 4);

%% 
%{ 
%
% --------------- Read In Quarterly Data for Nominal PCE ---------------- 
% xlsname = '/Users/mwatson/Dropbox/p_comp/revision/data/nom_pce_m59.xlsx';
xlsname = 'data/nom_pce_m59.xlsx';
% Read Data 
ndesc=1;
ncodes=1;
sheet='Sheet1';
[namevec,descmat,tcodemat,datevec,datamat] = readxls(xlsname,sheet,ns,dnobs_m,ndesc,ncodes);
tcodemat = tcodemat(1,:);
labelvec=descmat(:,1);
% Convert Namestrings to upper case 
namevec = upper(namevec);
% Eliminate any leading or trailing blanks 
namevec=strtrim(namevec);
labelvec=strtrim(labelvec);
% Replace missing values with NaN
isel = datamat == miss_code;
datamat(isel) = NaN;
pce_nom = datamat;
namevec_nom = namevec;
labelvec_nom = labelvec;
tcodemat_nom = tcodemat;

% ------------- Read in Price Deflators and compute inflation Rates ---------
% xlsname = '/Users/mwatson/Dropbox/p_comp/revision/data/p_pce_m59.xlsx';
xlsname = 'data/p_pce_m59.xlsx';
% Read Data 
ndesc=2;
ncodes=2;
sheet='Sheet1';
[namevec,descmat,codemat,datevec,datamat] = readxls(xlsname,sheet,ns,dnobs_m,ndesc,ncodes);
tcodemat = codemat(1,:);
ordermat = codemat(2,:)';
labelvec=descmat(:,1);
labelvec_short = descmat(:,2);
% Convert Namestrings to upper case 
namevec = upper(namevec);
% Eliminate any leading or trailing blanks 
namevec=strtrim(namevec);
labelvec=strtrim(labelvec);
labelvec_short=strtrim(labelvec_short);
% Replace missing values with NaN
isel = datamat == miss_code;
datamat(isel) = NaN;
pce_p = datamat;
namevec_p = namevec;
labelvec_p = labelvec;
labelvec_short_p = labelvec_short;
tcodemat_p = tcodemat;

% Compute Total and Values for Included Series
isel = tcodemat_p==2;
p_agg = pce_p(:,isel==1);         % Aggregrate Price Indenx  
isel = tcodemat==2;
nom_agg = pce_nom(:,isel==1);     % Aggregate Nominal PCE
isel = tcodemat_p==5;
p_agg_xfe = pce_p(:,isel==1);     % Aggregrate Price Index XFE  
isel = tcodemat==5;
nom_agg_xfe = pce_nom(:,isel==1);  % Aggregate Nominal PCE XFE

% Compute Energy Services Component of Household Expenditures
isel = tcodemat==3;
nom_eg = pce_nom(:,isel==1);   % Gasoline and other energy goods
isel = tcodemat==4;
nom_hu = pce_nom(:,isel==1);    % Housing and utilities
isel = tcodemat==6;
nom_egs = pce_nom(:,isel==1);   % Energy Goods and Services
isel = tcodemat_p==3;
p_eg = pce_p(:,isel==1);        % Gasoline and other energy goods
isel = tcodemat_p==4;
p_hu = pce_p(:,isel==1);        % Housing and utilities
isel = tcodemat_p==6;
p_egs = pce_p(:,isel==1);       % Energy Goods and Services

nom_hu_es = nom_egs - nom_eg;    % Energy services component of housing and utilities
nom_hu_xes = nom_hu - nom_hu_es; % Housing and Utilities excluding energy services

% Compute some inflation rates
infl_eg = NaN*zeros(dnobs_m,1);
infl_eg(2:end) = log(p_eg(2:end)./p_eg(1:end-1));
infl_egs = NaN*zeros(dnobs_m,1);
infl_egs(2:end) = log(p_egs(2:end)./p_egs(1:end-1));
s_eg_egs = nom_eg./nom_egs;   % Share of energy goods in energy goods and services 
infl_hu_es = (infl_egs - (s_eg_egs.*infl_eg))./(1-s_eg_egs);
p_hu_es = NaN*zeros(dnobs_m,1);
p_hu_es(1) = 100;
for i = 2:dnobs_m;
    p_hu_es(i) = p_hu_es(i-1)*exp(infl_hu_es(i));
end;
infl_hu = NaN*zeros(dnobs_m,1);
infl_hu(2:end) = log(p_hu(2:end)./p_hu(1:end-1));
s_es_hu = nom_hu_es./nom_hu;  % Share of energy in housing and utilities
infl_hu_xes = (infl_hu - (s_es_hu.*infl_hu_es))./(1-s_es_hu);
p_hu_xes = NaN*ones(dnobs_m,1);
p_hu_xes(1) = 100;
for i = 2:dnobs_m;
    p_hu_xes(i) = p_hu_xes(i-1)*exp(infl_hu_xes(i));
end;

% Add New Series to various matrices
pce_nom = [pce_nom nom_hu_es nom_hu_xes];
tcodemat = [tcodemat 1 1];
namevec = [namevec;'nom_hu_es';'nom_hu_xes'];
labelvec = [labelvec;'Energy Services in Housing and Utilities';'Housing and Utilities excl Energy Services'];

pce_p = [pce_p [p_hu_es p_hu_xes]];
tcodemat_p = [tcodemat_p 1 1];
ordermat_p = [ordermat;9.2;9.1];
namevec_p = [namevec_p;'p_hu_es';'p_hu_xes'];
labelvec_p = [labelvec_p;'Gas & electric untilies';'Housing excluding gas & electric utilities'];
labelvec_short_p = [labelvec_short_p;'Gas & electric untilies';'Housing excl. gas & elec. util.'];

isel = (tcodemat_p==1)+(tcodemat_p==3);
p_disagg = pce_p(:,isel==1);
ordervec = ordermat_p(isel'==1);
nom_disagg = pce_nom(:,isel==1);
namevec_disagg = namevec(isel==1);
labelvec_disagg = labelvec_p(isel==1);
labelvec_short_disagg = labelvec_short_p(isel==1);
n_incl = size(p_disagg,2);

% ... Data Check ... 
 % Add nominal disaggegrated PCE and compare to aggregate
 tmp = sum(nom_disagg,2);
 tmp_dif = (tmp - nom_agg)./nom_agg; 
 if max(abs(tmp_dif)) > .0001;
     error('nominal dissaggregates do not add to nominal aggregate');
 end;
 

%}


%% Load GT data

% Loads GT data = p_gt 
load("data/gt_data.mat") 

% Load raw data series
p_agg = p_gt(1:end-2, 1); 
p_agg_xfe = p_gt(1:end-2, 2); 
p_agg_xe = p_gt(1:end-2, 3); 

 % Compute Various Inflation Measures for Aggregate and Dissaggregate PCE
 % -- Monthly data
 p_agg_m = p_agg; % Monthly aggregate index
 p_agg_xfe_m = p_agg_xfe; % Monthly ex food & energy 
 p_agg_xe_m = p_agg_xe; % Monthly ex energy 

 % For disaggregated data
 % p_disagg_m = p_disagg; % Monthly disaggregated data
 

 %{
 nom_disagg_m = nom_disagg;
 nom_agg_m = nom_agg;
 nom_agg_xfe_m = nom_agg_xfe;
 pce_nom_m = pce_nom;
 nom_egs_m = nom_egs;
 p_egs_m = p_egs;
 %}


 % -- Quarterly data
 p_agg_q = mtoq(p_agg,calds_m,calds_q,mtoq_agg);
 p_agg_xfe_q = mtoq(p_agg_xfe,calds_m,calds_q,mtoq_agg);
 p_agg_xe_q = mtoq(p_agg_xe,calds_m,calds_q,mtoq_agg);
 
 % p_disagg_q = mtoq(p_disagg,calds_m,calds_q,mtoq_agg);
 
 %{
 nom_disagg_q = mtoq(nom_disagg,calds_m,calds_q,mtoq_agg);
 nom_agg_q = mtoq(nom_agg,calds_m,calds_q,mtoq_agg);
 nom_agg_xfe_q = mtoq(nom_agg_xfe,calds_m,calds_q,mtoq_agg);
 pce_nom_q = mtoq(pce_nom,calds_m,calds_q,mtoq_agg);
 nom_egs_q = mtoq(nom_egs,calds_m,calds_q,mtoq_agg);
 p_egs_q = mtoq(p_egs,calds_m,calds_q,mtoq_agg);
 %}

 % --- Monthly Inflation ----
 paar = 1200;
 dp_agg_m = paar*dif(log(p_agg_m),1);
 dp_agg_xfe_m = paar*dif(log(p_agg_xfe_m),1);
 dp_agg_xe_m = paar*dif(log(p_agg_xe_m),1);
 % dp_disagg_m = paar*dif(log(p_disagg_m),1);
 
 % --- Quarterly Inflation ----
 paar = 400;
 dp_agg_q = paar*dif(log(p_agg_q),1);
 dp_agg_xfe_q = paar*dif(log(p_agg_xfe_q),1);
 dp_agg_xe_q = paar*dif(log(p_agg_xe_q),1);
 % dp_disagg_q = paar*dif(log(p_disagg_q),1);

%%

%{
 % Compute dp_agg_xe:  aggregate inflation excluding energy
 % -- Monthly
 s_egs_agg_m = nom_egs_m./nom_agg_m;            % Share of energy goods and services in total
 dp_egs_m = 1200*dif(log(p_egs_m),1);           % Energy goods and services Inflation
 dp_agg_xe_m = (dp_agg_m - s_egs_agg_m.*dp_egs_m)./(1-s_egs_agg_m);
 tmp = [0;dp_agg_xe_m(2:end)/1200];
 tmp = cumsum(tmp);
 p_agg_xe_m = exp(tmp);
 
 % -- Quarterly
 p_agg_xe_q = mtoq(p_agg_xe_m,calds_m,calds_q,mtoq_agg);
 dp_agg_xe_q = 400*dif(log(p_agg_xe_q),1);
 s_egs_agg_q = nom_egs_q./nom_agg_q;            % Share of energy goods and services in total
 dp_egs_q = 400*dif(log(p_egs_q),1);            % Energy goods and services Inflation
 % dp_agg_xe_q = (dp_agg_q - s_egs_agg_q.*dp_egs_q)./(1-s_egs_agg_q); 



%  % ---- Plot Series
%  figure;
%  subplot(1,2,1);
%  plot(calvec_m,dp_agg_m);
%  title('Headline Inflation');
%  subplot(1,2,2);
%  plot(calvec_q,dp_agg_q);
%  title('Headline Inflation');
%  waitforbuttonpress;
%  
%  figure;
%  subplot(1,2,1);
%  plot(calvec_m,dp_agg_xfe_m);
%  title('Core-XFE Inflation');
%  subplot(1,2,2);
%  plot(calvec_q,dp_agg_xfe_q);
%  title('Core-XFE Inflation');
%  waitforbuttonpress;
%  
%  figure;
%  subplot(1,2,1);
%  plot(calvec_m,dp_agg_xe_m);
%  title('Core-XE Inflation');
%  subplot(1,2,2);
%  plot(calvec_q,dp_agg_xe_q);
%  title('Core-XE Inflation');
%  waitforbuttonpress;
%  
%  for i = 1:n_incl;
%      str = labelvec_disagg(i);
%      figure;
%      subplot(1,2,1);
%      plot(calvec_m,dp_disagg_m(:,i));
%      title(str);
%      subplot(1,2,2);
%      plot(calvec_q,dp_disagg_q(:,i));
%      title(str);
%      waitforbuttonpress;
%  end;
 
 %  Compute Shares of nominal PCE
 %  --- Monthly
 tmp = sum(nom_disagg_m,2);
 share_disagg_m = nom_disagg_m./repmat(tmp,1,n_incl);
 share_avg_m = NaN*zeros(dnobs_m,n_incl);
 share_avg_m(2:end,:) = (nom_disagg_m(2:end,:)+nom_disagg_m(1:end-1,:))./(repmat(tmp(2:end),1,n_incl)+repmat(tmp(1:end-1),1,n_incl));
 
 %  --- Quarterly
 tmp = sum(nom_disagg_q,2);
 share_disagg_q = nom_disagg_q./repmat(tmp,1,n_incl);
 share_avg_q = NaN*zeros(dnobs_q,n_incl);
 share_avg_q(2:end,:) = (nom_disagg_q(2:end,:)+nom_disagg_q(1:end-1,:))./(repmat(tmp(2:end),1,n_incl)+repmat(tmp(1:end-1),1,n_incl));
 
 
 % -- Compute Shares for XFE Components
 i_xfe = [5 7 16];
 i_xe = [7 16];
 share_avg_xfe_m = share_avg_m;
 share_avg_xfe_m(:,i_xfe) = 0;
 share_avg_xfe_m = share_avg_xfe_m./repmat(sum(share_avg_xfe_m,2),1,size(share_avg_m,2));
 share_avg_xe_m = share_avg_m;
 share_avg_xe_m(:,i_xe) = 0;
 share_avg_xe_m = share_avg_xe_m./repmat(sum(share_avg_xe_m,2),1,size(share_avg_m,2));
 share_avg_xfe_q = share_avg_q;
 share_avg_xfe_q(:,i_xfe) = 0;
 share_avg_xfe_q = share_avg_xfe_q./repmat(sum(share_avg_xfe_q,2),1,size(share_avg_q,2));
 share_avg_xe_q = share_avg_q;
 share_avg_xe_q(:,i_xe) = 0;
 share_avg_xe_q = share_avg_xe_q./repmat(sum(share_avg_xe_q,2),1,size(share_avg_q,2));

%  % Plot Shares
%  for i = 1:n_incl;
%      str = labelvec_disagg(i);
%      figure;
%      subplot(1,2,1)
%       plot(calvec_m,[share_disagg_m(:,i) share_avg_m(:,i)]); 
%       title(str);
%      subplot(1,2,2);
%       plot(calvec_q,[share_disagg_q(:,i) share_avg_q(:,i)]); 
%       title(str);
%      waitforbuttonpress;
%  end;
 
% Compute Approximations to Aggregate inflation
 % -- Monthly -- 
 tmp = dp_disagg_m.*share_avg_m;
 dp_agg_approx1_m = sum(tmp,2);
 weight_avg_m = mean(share_disagg_m,1)';
 dp_agg_approx2_m = dp_disagg_m*weight_avg_m;
 dif1_m = dp_agg_m - dp_agg_approx1_m;
 dif2_m = dp_agg_m - dp_agg_approx2_m;
 
 % -- Quarterly -- 
 tmp = dp_disagg_q.*share_avg_q;
 dp_agg_approx1_q = sum(tmp,2);
 weight_avg_q = mean(share_disagg_q,1)';
 dp_agg_approx2_q = dp_disagg_q*weight_avg_q;
 dif1_q = dp_agg_q - dp_agg_approx1_q;
 dif2_q = dp_agg_q - dp_agg_approx2_q;
 

%  figure;
%  subplot(1,2,1);
%   plot(calvec_m,[dp_agg_m dp_agg_approx1_m dif1_m]);
%  subplot(1,2,2);
%   plot(calvec_q,[dp_agg_q dp_agg_approx1_q dif1_q]);
%  waitforbuttonpress;
%  figure;
%  subplot(1,2,1);
%    plot(calvec_m,[dp_agg_m dp_agg_approx2_m dif2_m]);
%  subplot(1,2,2);
%    plot(calvec_q,[dp_agg_q dp_agg_approx2_q dif2_q]);
%  waitforbuttonpress;


% Reorder components by ordervec;
%[tmp,i_order] = sort(weight_avg,'descend');
[tmp,i_order] = sort(ordervec,'ascend');
nom_disagg_m = nom_disagg_m(:,i_order);
dp_disagg_m = dp_disagg_m(:,i_order);
share_avg_m = share_avg_m(:,i_order);
share_avg_xe_m = share_avg_xe_m(:,i_order);
share_avg_xfe_m = share_avg_xfe_m(:,i_order);
weight_avg_m = weight_avg_m(i_order);
nom_disagg_q = nom_disagg_q(:,i_order);
dp_disagg_q = dp_disagg_q(:,i_order);
share_avg_q = share_avg_q(:,i_order);
share_avg_xe_q = share_avg_xe_q(:,i_order);
share_avg_xfe_q = share_avg_xfe_q(:,i_order);
weight_avg_q = weight_avg_q(i_order);
labelvec_disagg = labelvec_disagg(i_order);
labelvec_short_disagg = labelvec_short_disagg(i_order);
namevec_disagg = namevec_disagg(i_order);

% Compute a new dissagregated version of inflation using cor_xfe, food, and
% energy
str='DFXARG';    % food
str=upper(str);
j = colnumber(str,namevec_disagg);
% -- Monthly -- 
dp_food_m = dp_disagg_m(:,j);
nom_food_m = nom_disagg_m(:,j);
isel = tcodemat==5;
nom_agg_xfe_m = pce_nom_m(:,isel==1);     % Nominal PCE 
nom_disagg_3comp_m = [nom_agg_xfe_m nom_food_m nom_egs_m];
tmp = sum(nom_disagg_3comp_m,2);
share_disagg_3comp_m = nom_disagg_3comp_m./repmat(tmp,1,3);
share_avg_3comp_m = NaN*zeros(dnobs_m,3);
share_avg_3comp_m(2:end,:) = (nom_disagg_3comp_m(2:end,:)+nom_disagg_3comp_m(1:end-1,:))./(repmat(tmp(2:end),1,3)+repmat(tmp(1:end-1),1,3));
weight_avg_3comp_m = mean(share_disagg_3comp_m,1)';
dp_disagg_3comp_m = [dp_agg_xfe_m dp_food_m dp_egs_m];
% -- Quarterly -- 
dp_food_q = dp_disagg_q(:,j);
nom_food_q = nom_disagg_q(:,j);
isel = tcodemat==5;
nom_agg_xfe_q = pce_nom_q(:,isel==1);     % Nominal PCE 
nom_disagg_3comp_q = [nom_agg_xfe_q nom_food_q nom_egs_q];
tmp = sum(nom_disagg_3comp_q,2);
share_disagg_3comp_q = nom_disagg_3comp_q./repmat(tmp,1,3);
share_avg_3comp_q = NaN*zeros(dnobs_q,3);
share_avg_3comp_q(2:end,:) = (nom_disagg_3comp_q(2:end,:)+nom_disagg_3comp_q(1:end-1,:))./(repmat(tmp(2:end),1,3)+repmat(tmp(1:end-1),1,3));
weight_avg_3comp_q = mean(share_disagg_3comp_q,1)';
dp_disagg_3comp_q = [dp_agg_xfe_q dp_food_q dp_egs_q];

share_avg_3comp_xfe_m = share_avg_3comp_m;
share_avg_3comp_xfe_m(:,[2 3]) = 0;
share_avg_3comp_xfe_m(:,1) = 1;
share_avg_3comp_xfe_q = share_avg_3comp_q;
share_avg_3comp_xfe_q(:,[2 3]) = 0;
share_avg_3comp_xfe_q(:,1) = 1;
share_avg_3comp_xe_m = share_avg_3comp_m;
share_avg_3comp_xe_m(:,3) = 0;
share_avg_3comp_xe_m = share_avg_3comp_xe_m./repmat(sum(share_avg_3comp_xe_m,2),1,size(share_avg_3comp_m,2));
share_avg_3comp_xe_q = share_avg_3comp_q;
share_avg_3comp_xe_q(:,3) = 0;
share_avg_3comp_xe_q = share_avg_3comp_xe_q./repmat(sum(share_avg_3comp_xe_q,2),1,size(share_avg_3comp_q,2));

labelvec_disagg_3comp = ['CoreXFE';'Food   ';'Energy '];
namevec_disagg_3comp = ['CoreXFE';'Food   ';'Energy '];
n_incl_3comp = 3;
%}

%{
 % Save Variable Series 
  slist = {...
           'nom_agg_m' ... 
           'nom_disagg_m' ...   
           'dp_agg_m' ...
           'dp_agg_xfe_m' ...
           'dp_agg_xe_m' ...
           'dp_disagg_m' ...
           'dp_agg_approx1_m' ...
           'dp_agg_approx2_m' ...
           'share_avg_m' ...
           'weight_avg_m' ...
           'dp_disagg_3comp_m' ...
           'share_avg_3comp_m' ...
           'weight_avg_3comp_m' ...
           'calvec_m' ...
           'calds_m' ...
           'dnobs_m' ...
           'nom_agg_q' ... 
           'nom_disagg_q' ...   
           'dp_agg_q' ...
           'dp_agg_xfe_q' ...
           'dp_agg_xe_q' ...
           'dp_disagg_q' ...
           'dp_agg_approx1_q' ...
           'dp_agg_approx2_q' ...
           'share_avg_q' ...
           'weight_avg_q' ...
           'dp_disagg_3comp_q' ...
           'share_avg_3comp_q' ...
           'weight_avg_3comp_q' ...
           'calvec_q' ...
           'calds_q' ...
           'dnobs_q' ...
           'n_incl' ...
           'labelvec_disagg' ...
           'labelvec_short_disagg' ...
           'namevec_disagg' ...
           'labelvec_disagg_3comp' ...
           'namevec_disagg_3comp' ...
           'n_incl_3comp' ...
           'share_avg_xfe_m' ...
           'share_avg_xe_m' ...
           'share_avg_xfe_q' ...
           'share_avg_xe_q' ...
           'share_avg_3comp_xfe_m' ...
           'share_avg_3comp_xe_m' ...
           'share_avg_3comp_xfe_q' ...
           'share_avg_3comp_xe_q' ...
          }';
  str_tmp = [matdir 'slist'];
  save(str_tmp,'slist');
  for iseries = 1:size(slist,1);
     ustr = char(slist(iseries));
     str_tmp = [matdir ustr];
     str_tmp1 = 'save(str_tmp,''';
     str_tmp2 = ''');';
     eval([str_tmp1 ustr str_tmp2]);
  end;

  
% Load Variables and Give Standard Names 
str_tmp = [matdir 'slist'];
load(str_tmp); 
for iseries = 1:size(slist,1);
     ustr = char(slist(iseries));
     str_tmp = [matdir ustr];
     load(str_tmp);
end; 

%}

end;  %Load Data end;