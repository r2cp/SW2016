% Table 3
% November 4, 2015
%


clear all;
small = 1.0e-10;
big = 1.0e+6;                  
rng(63761);
  
  % -- File Directories  
  outdir = 'matlab/out/';
  figdir = 'matlab/fig/';
  matdir = 'matlab/mat/gtM/';
  matdir_nojumps = 'matlab/mat/gtM-nojumps/';
  matdir_constvol = 'matlab/mat/gtM-const_vol_e/';
  
  % -- Read in Data --- 
  load_data = 1;  % 1 if reloading data from Excel, etc 
  mtoq_agg = 0;   % Temporal aggregation indicator of monthly to quarterly data
  nper = 12;
  pcomp_data_calendar_m_and_q_gt;

  % Data Series Used
  dp_agg = dp_agg_m;
  dp_agg_xfe = dp_agg_xfe_m;
  dp_agg_xe = dp_agg_xe_m;
%   dp_disagg = dp_disagg_q;
%   share_avg = share_avg_q;
  
  calvec = calvec_m;
  dnobs = dnobs_m;
  calds = calds_m;
  
  label_suffix = '_gtM';
  
  % Dates
  first_date = [2001 1];
  last_date = calds(end,:);
  ismpl_60_end = smpl(calvec,first_date,last_date,nper);
  cal_60_end = calvec(ismpl_60_end==1);
  calds_60_end = calds(ismpl_60_end==1,:);
  dnobs_60_end = size(cal_60_end,1);
  
  % First POOS Date
  first_date_poos = [2015 1];
  tmp = calds_60_end == repmat(first_date_poos,dnobs_60_end,1);
  [a,i_first_poos] = max(sum(tmp,2));  % Index in calvec_ismpl of first poos period
  
%   % Date for 2008:Q4 -- skip
%   if nper == 4;
%    date_skip = [2008 4];
%    tmp = calds_60_end == repmat(date_skip,dnobs_60_end,1);
%    [a,i_skip] = max(sum(tmp,2));  % Index in calvec_ismpl of first poos period
%   else; 
%    date_skip = [2008 10];
%    tmp = calds_60_end == repmat(date_skip,dnobs_60_end,1);
%    [a,ii] = max(sum(tmp,2));  % Index in calvec_ismpl of first poos period
%    i_skip = [ii ii+1 ii+2]';
%   end;
  
  % Create output file
  outfile_name = [outdir 'fcst_rslt_poos' label_suffix '.txt'];
  fileID = fopen(outfile_name,'w');
  
  %% Used for lagged models (contemporaneous values as forecast) 
  % Slicing for the periods with available values
  dp_agg_60_end = dp_agg(ismpl_60_end==1,:);
  dp_agg_xe_60_end = dp_agg_xe(ismpl_60_end==1,:);
  dp_agg_xfe_60_end = dp_agg_xfe(ismpl_60_end==1,:);
%   dp_disagg_60_end = dp_disagg(ismpl_60_end==1,:);
%   share_avg_60_end = share_avg(ismpl_60_end==1,:);
  

%% Atkeson-Ohanian Random walk models of inflation 
 
  dpa_agg = zeros(dnobs,1);
  dpa_agg_xe = zeros(dnobs,1);
  dpa_agg_xfe = zeros(dnobs,1);
%   dpa_disagg =  zeros(dnobs,n_incl);
  for i = 1:nper;
    dpa_agg(nper:end) = dpa_agg(nper:end)+dp_agg(i:end-nper+i);
    dpa_agg_xe(nper:end) = dpa_agg_xe(nper:end)+dp_agg_xe(i:end-nper+i);
    dpa_agg_xfe(nper:end) = dpa_agg_xfe(nper:end)+dp_agg_xfe(i:end-nper+i);
%     dpa_disagg(nper:end,:) = dpa_disagg(nper:end,:)+dp_disagg(i:end-nper+i,:);
  end;
  dpa_agg = dpa_agg/nper;
  dpa_agg_xe = dpa_agg_xe/nper;
  dpa_agg_xfe = dpa_agg_xfe/nper;
%   dpa_disagg = dpa_disagg/nper;

  % Slicing for the relevant sample 
  dpa_agg_60_end = dpa_agg(ismpl_60_end==1);
  dpa_agg_xe_60_end = dpa_agg_xe(ismpl_60_end==1);
  dpa_agg_xfe_60_end = dpa_agg_xfe(ismpl_60_end==1);
%   dpa_disagg_60_end = dpa_disagg(ismpl_60_end==1,:);
  
  %% -- Univariate Filtered and Smoothed Values

  % Removed the _poos for now
  % ulabel = ['dp_agg_poos' label_suffix];

  % Baseline results for monthly model 
  ulabel = ['dp_agg' label_suffix];
  str_tmp = [matdir 'tau_mean_pct' ulabel]; load(str_tmp);
  tau_f_h = tau_mean_pct(:,1);
  
  ulabel = ['dp_xe' label_suffix];
  str_tmp = [matdir 'tau_mean_pct' ulabel]; load(str_tmp);
  tau_f_xe = tau_mean_pct(:,1);
  
  ulabel = ['dp_xfe' label_suffix];
  str_tmp = [matdir 'tau_mean_pct' ulabel]; load(str_tmp);
  tau_f_xfe = tau_mean_pct(:,1);


  % No outlier-adjustment results
  ulabel = ['dp_agg' label_suffix];
  str_tmp = [matdir_nojumps 'tau_mean_pct' ulabel]; load(str_tmp);
  tau_f_h_nojumps = tau_mean_pct(:,1);
  
  ulabel = ['dp_xe' label_suffix];
  str_tmp = [matdir_nojumps 'tau_mean_pct' ulabel]; load(str_tmp);
  tau_f_xe_nojumps = tau_mean_pct(:,1);
  
  ulabel = ['dp_xfe' label_suffix];
  str_tmp = [matdir_nojumps 'tau_mean_pct' ulabel]; load(str_tmp);
  tau_f_xfe_nojumps = tau_mean_pct(:,1);

  % Results from constant volatility of transitory component 
  ulabel = ['dp_agg' label_suffix];
  str_tmp = [matdir_constvol 'tau_mean_pct' ulabel]; load(str_tmp);
  tau_f_h_constvol = tau_mean_pct(:,1);
  
  ulabel = ['dp_xe' label_suffix];
  str_tmp = [matdir_constvol 'tau_mean_pct' ulabel]; load(str_tmp);
  tau_f_xe_constvol = tau_mean_pct(:,1);
  
  ulabel = ['dp_xfe' label_suffix];
  str_tmp = [matdir_constvol 'tau_mean_pct' ulabel]; load(str_tmp);
  tau_f_xfe_constvol = tau_mean_pct(:,1);



%   % -- Some results for the 17-component model
%   mlabel = ['_mucsv_tvp_17c' label_suffix];
%   str_tmp = [matdir 'agg_tau_mean_pct_poos' mlabel]; load(str_tmp);
%   tau_f_17c = agg_tau_mean_pct(:,1);
%   tau_f_mean_pct_17c = agg_tau_mean_pct;
%   mlabel = ['_mucsv_constAlpha_17c' label_suffix];
%   str_tmp = [matdir 'agg_tau_mean_pct_poos' mlabel]; load(str_tmp);
%   tau_f_constAlpha_17c = agg_tau_mean_pct(:,1);
%   tau_f_constAlpha_mean_pct_17c = agg_tau_mean_pct;
% %  tau_f_17c = tau_f_h;
  
%   % -- Some results for the 3-component model
%   mlabel = ['_mucsv_tvp_3c' label_suffix];
%   str_tmp = [matdir 'agg_tau_mean_pct_poos' mlabel]; load(str_tmp);
%   tau_f_3c = agg_tau_mean_pct(:,1);
%   tau_f_mean_pct_3c = agg_tau_mean_pct;
%   mlabel = ['_mucsv_constAlpha_3c' label_suffix];
%   str_tmp = [matdir 'agg_tau_mean_pct_poos' mlabel]; load(str_tmp);
%   tau_f_constAlpha_3c = agg_tau_mean_pct(:,1);
%   tau_f_constAlpha_mean_pct_3c = agg_tau_mean_pct;
% %  tau_f_3c = tau_f_h;
  
  %% Compute target from average values of CPI-all 

  nhorizon_years = [1 2 3];

  for ih = 1: size(nhorizon_years,2)
   
    nhorizon_avg = nhorizon_years(ih) * nper;
    % Compute Target series to forecast
    dp_avg_target = NaN(dnobs_60_end,1);
    dp_avg_target_x = NaN(dnobs_60_end,1); % Used for panel (B), where SW exclude a specific period
  
    for t = i_first_poos:dnobs_60_end-nhorizon_avg
      % Target is average inflation from t+1 to t+h, where h=12,24,36 for
      % monthly observations 
      dp_avg_target(t) = mean(dp_agg_60_end(t+1:t+nhorizon_avg));
    end
     
     % For the special 2008 skip in SW: computes a 2nd column for the
     % errors
%      dp_x = dp_agg_60_end;
%      dp_x(i_skip) = NaN;
%      for t = i_first_poos:dnobs_60_end-nhorizon_avg;
%        tmp = dp_x(t+1:t+nhorizon_avg);
%        dp_avg_target_x(t) = mean(tmp(isnan(tmp)==0));
%      end;
  
    % Compute Forecast errors
    % target = [dp_avg_target dp_avg_target_x];
    target = [dp_avg_target dp_avg_target]; % Same column twice
    
    % Baseline model results
    e_h = repmat(tau_f_h,1,2)- target;
    e_xe = repmat(tau_f_xe,1,2)- target;
    e_xfe = repmat(tau_f_xfe,1,2)- target;

    % No outlier-adjustment results
    e_h_nojumps = repmat(tau_f_h_nojumps,1,2) - target;
    e_xe_nojumps = repmat(tau_f_xe_nojumps,1,2) - target;
    e_xfe_nojumps = repmat(tau_f_xfe_nojumps,1,2) - target;

    % No outlier-adjustment results
    e_h_constvol = repmat(tau_f_h_constvol,1,2) - target;
    e_xe_constvol = repmat(tau_f_xe_constvol,1,2) - target;
    e_xfe_constvol = repmat(tau_f_xfe_constvol,1,2) - target;

%      e_17c = repmat(tau_f_17c,1,2)- target;
%      e_3c = repmat(tau_f_3c,1,2)- target;
%      e_constAlpha_17c = repmat(tau_f_constAlpha_17c,1,2)- target;
%      e_constAlpha_3c = repmat(tau_f_constAlpha_3c,1,2)- target;

    % Lagged models
     e_l = repmat(dp_agg_60_end,1,2)- target;
     e_l_xe = repmat(dp_agg_xe_60_end,1,2)- target;
     e_l_xfe = repmat(dp_agg_xfe_60_end,1,2)- target;

     % Atkeson-Ohanian random walk models
     e_a = repmat(dpa_agg_60_end,1,2)- target;
     e_a_xe = repmat(dpa_agg_xe_60_end,1,2)- target;
     e_a_xfe = repmat(dpa_agg_xfe_60_end,1,2)- target;
  
     % Compute MSFE = Mean Squared Forecast Error  
  
     mse = NaN(15,2);
     se_mse = NaN(15,2);
     n_mse = NaN(15,2);
     dif_mse = NaN(15,2);
     se_dif_mse = NaN(15,2);
  
     nma = nhorizon_avg+2;
     for j = 1:2;

       % Columns of error series from all the models 
       % e matrix with errors for each column 
       % e = [e_17c(:,j) e_3c(:,j) e_constAlpha_17c(:,j) e_constAlpha_3c(:,j) e_h(:,j) e_xe(:,j) e_xfe(:,j) e_l(:,j) e_l_xe(:,j) e_l_xfe(:,j) e_a(:,j) e_a_xe(:,j) e_a_xfe(:,j)];
       
       e = [... 
            ... % Baseline results
            e_h(:,j) e_xe(:,j) e_xfe(:,j) ... 
            ... % No outlier-adjustment results
            e_h_nojumps(:,j) e_xe_nojumps(:,j) e_xfe_nojumps(:,j) ... 
            ... % Constant volatility results
            e_h_constvol(:,j) e_xe_constvol(:,j) e_xfe_constvol(:,j) ... 
            ... % Lagged values 
            e_l(:,j) e_l_xe(:,j) e_l_xfe(:,j) ... 
            ... % RW models 
            e_a(:,j) e_a_xe(:,j) e_a_xfe(:,j)];
       
       % For each of the error columns: 
       for i = 1:size(e,2);
        tmp = e(:,i);
        tmp1 = e(:,1); % 1st column is used to measure the differences
        
        % Get non NaN positions 
        tmp = tmp(isnan(tmp) == 0);
        tmp1 = tmp1(isnan(tmp1) == 0);
        n_mse(i,j) = size(tmp,1);
        
        % Compute MSE and its standard error
        tmp = tmp.^2;
        tmp1 = tmp1.^2;   
        [b,vb] = hacm(tmp,nma,1);
        mse(i,j) = b;
        % mse(i,j) = sqrt(mean(tmp)); % RMSFE 
        se_mse(i,j) = sqrt(vb);
        % se_mse(i,j) = sqrt(sqrt(vb)); % std error of RMSFE

        % Compute difference and its standard error
        tmpdif = tmp-tmp1;
        [b,vb] = hacm(tmpdif,nma,1);
        dif_mse(i,j) = b;
        % dif_mse(i,j) = sqrt(mean(tmpdif));
        se_dif_mse(i,j) = sqrt(vb);
       end
     end
     
     % If for the first horizon, define the output matrices
     if ih == 1
         mse_mat = NaN(size(mse,1),size(mse,2),size(nhorizon_years,2));
         se_mse_mat = NaN(size(mse,1),size(mse,2),size(nhorizon_years,2));
         dif_mse_mat = NaN(size(mse,1),size(mse,2),size(nhorizon_years,2));
         se_dif_mse_mat = NaN(size(mse,1),size(mse,2),size(nhorizon_years,2));
     end
     mse_mat(:,:,ih) = mse;
     se_mse_mat(:,:,ih) = se_mse;
     dif_mse_mat(:,:,ih) = dif_mse;
     se_dif_mse_mat(:,:,ih) = se_dif_mse;
         
  end
  
  %% Output file with forecast evaluation 
  fprintf(fileID,'All Observations \n'); 
  % In this loop i = model used 
  for i = 1:size(e,2);
    for ih = 1:size(nhorizon_years,2);
      fprintf(fileID,'%4.2f ',mse_mat(i,1,ih));
      fprintf(fileID,'(%4.2f),',se_mse_mat(i,1,ih));
      fprintf(fileID,'%4.2f ',dif_mse_mat(i,1,ih));
      
      if ih < size(nhorizon_years,2)
       fprintf(fileID,'(%4.2f),',se_dif_mse_mat(i,1,ih));
      else
        fprintf(fileID,'(%4.2f) \n',se_dif_mse_mat(i,1,ih));
      end
    end  
  end
  
%   fprintf(fileID,'\n\n Excluding 2008:Q4 \n');  
%   for i = 1:size(e,2);
%     for ih = 1:size(nhorizon_years,2);
%       fprintf(fileID,'%4.2f ',mse_mat(i,2,ih));
%       fprintf(fileID,'(%4.2f),',se_mse_mat(i,2,ih));
%       fprintf(fileID,'%4.2f ',dif_mse_mat(i,2,ih));
%       if ih < size(nhorizon_years,2);
%        fprintf(fileID,'(%4.2f),',se_dif_mse_mat(i,2,ih));
%       else;
%         fprintf(fileID,'(%4.2f) \n',se_dif_mse_mat(i,2,ih));
%       end;
%     end;   
%   end;
  
  

 