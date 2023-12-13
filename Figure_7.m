% Figure 7 
% November 4, 2015
%


clear all;
small = 1.0e-10;
big = 1.0e+6;                  
rng(63761);
 
  % -- File Directories  
  outdir = '/Users/mwatson/Dropbox/p_comp/revision/matlab/out/';
  figdir = '/Users/mwatson/Dropbox/p_comp/revision/matlab/fig/';
  matdir = '/Users/mwatson/Dropbox/p_comp/revision/matlab/mat/';
  
    % -- Read in Data --- 
  load_data = 1;  % 1 if reloading data from Excel, etc 
  mtoq_agg = 3;   % Temporal aggregation indicator of monthly to quarterly data
  pcomp_data_calendar_m_and_q;
  % Data Series Used
  dp_agg = dp_agg_q;
  dp_agg_xfe = dp_agg_xfe_q;
  dp_agg_xe = dp_agg_xe_q;
  dp_disagg = dp_disagg_q;
  share_avg = share_avg_q;
  calvec = calvec_q;
  dnobs = dnobs_q;
  calds = calds_q;
  nper = 4;
  mlabel = '_mucsv_tvp_17c_Q3';
  
  % Dates
  first_date = [1960 1];
  last_date = calds(end,:);
  ismpl_60_end = smpl(calvec,first_date,last_date,4);
  cal_60_end = calvec(ismpl_60_end==1);
  dnobs_60_end = size(cal_60_end,1);
  
  dp_agg_60_end = dp_agg(ismpl_60_end==1,:);
  dp_agg_xe_60_end = dp_agg_xe(ismpl_60_end==1,:);
  dp_agg_xfe_60_end = dp_agg_xfe(ismpl_60_end==1,:);
  dp_disagg_60_end = dp_disagg(ismpl_60_end==1,:);
  share_avg_60_end = share_avg(ismpl_60_end==1,:);
  
 
  str_tmp = [matdir 'agg_tau_mean_sd' mlabel]; load(str_tmp);
  tau_s_17c = agg_tau_mean_sd(:,1);
  str_tmp = [matdir 'agg_tau_mean_pct' mlabel]; load(str_tmp);
  tau_s_mean_pct_17c = agg_tau_mean_pct;
  str_tmp = [matdir 'disagg_tau_mean' mlabel]; load(str_tmp);
  dis_tau_s_17c = disagg_tau_mean;


  % -- Mean Values of Parameters
  str_tmp = [matdir 'alpha_eps_mean' mlabel]; load(str_tmp);
  str_tmp = [matdir 'alpha_tau_mean' mlabel]; load(str_tmp);
  str_tmp = [matdir 'alpha_eps_mean_pct' mlabel]; load(str_tmp);
  str_tmp = [matdir 'alpha_tau_mean_pct' mlabel]; load(str_tmp);
  str_tmp = [matdir 'sigma_eps_common_mean' mlabel]; load(str_tmp);
  str_tmp = [matdir 'sigma_eps_common_total_mean' mlabel]; load(str_tmp);
  str_tmp = [matdir 'sigma_eps_unique_mean' mlabel]; load(str_tmp);
  str_tmp = [matdir 'sigma_eps_unique_total_mean' mlabel]; load(str_tmp); 
  str_tmp = [matdir 'sigma_dtau_common_mean' mlabel]; load(str_tmp); 
  str_tmp = [matdir 'sigma_dtau_unique_mean' mlabel]; load(str_tmp);
  str_tmp = [matdir 'var_eps_common_mean' mlabel]; load(str_tmp);
  str_tmp = [matdir 'var_eps_common_total_mean' mlabel]; load(str_tmp);
  str_tmp = [matdir 'var_eps_unique_mean' mlabel]; load(str_tmp);
  str_tmp = [matdir 'var_eps_unique_total_mean' mlabel]; load(str_tmp); 
  str_tmp = [matdir 'var_y_eps_unique_mean_pct' mlabel]; load(str_tmp); 
  str_tmp = [matdir 'var_dtau_common_mean' mlabel]; load(str_tmp); 
  str_tmp = [matdir 'var_dtau_unique_mean' mlabel]; load(str_tmp);
  alpha_eps_17c = alpha_eps_mean;
  alpha_tau_17c = alpha_tau_mean;
  alpha_eps_pct_17c = alpha_eps_mean_pct;
  alpha_tau_pct_17c = alpha_tau_mean_pct;
  sigma_eps_common_17c = sigma_eps_common_mean;
  sigma_eps_common_total_17c = sigma_eps_common_total_mean;
  sigma_eps_unique_17c = sigma_eps_unique_mean;
  sigma_eps_unique_total_17c = sigma_eps_unique_mean;
  sigma_dtau_common_17c = sigma_dtau_common_mean;
  sigma_dtau_unique_17c = sigma_dtau_unique_mean;
  
  var_eps_common_17c = var_eps_common_mean;
  var_eps_common_total_17c = var_eps_common_total_mean;
  var_eps_unique_17c = var_eps_unique_mean;
  var_eps_unique_total_17c = var_eps_unique_total_mean;
  var_eps_unique_pct_17c = var_y_eps_unique_mean_pct;
  var_dtau_common_17c = var_dtau_common_mean;
  var_dtau_unique_17c = var_dtau_unique_mean;
  

% --- Figure 7 ----  

figure;
j = 14;
sigma_eps_j = sqrt(reshape(var_eps_unique_pct_17c(:,:,j),dnobs_60_end,6));

subplot(4,2,1);
   plot(cal_60_end,dp_disagg_60_end(:,j),'- k','LineWidth',1);
   tmp = char(labelvec_short_disagg(j));
  titstr =['(1) ' tmp]; 
  title(titstr);
  ax = gca;
  ax.FontSize = 20;
  ax.XLim = [1955 2020];
subplot(4,2,2);
  plot(cal_60_end,sigma_eps_j(:,1),'- k','LineWidth',3);
   hold on;
     plot(cal_60_end,sigma_eps_j(:,3),'-- k','LineWidth',3);
     plot(cal_60_end,sigma_eps_j(:,5),'-- k','LineWidth',3);
   hold off;
  title('\sigma_{\epsilon,{\iti}}');
  ax = gca;
  ax.FontSize = 20;
  ax.XLim = [1955 2020];
  
j = 5;
sigma_eps_j = sqrt(reshape(var_eps_unique_pct_17c(:,:,j),dnobs_60_end,6));

subplot(4,2,3);
   plot(cal_60_end,dp_disagg_60_end(:,j),'- k','LineWidth',1);
   tmp = char(labelvec_short_disagg(j));
  titstr =['(2) ' tmp]; 
  title(titstr);
  ax = gca;
  ax.FontSize = 20;
  ax.XLim = [1955 2020];
subplot(4,2,4);
  plot(cal_60_end,sigma_eps_j(:,1),'- k','LineWidth',3);
   hold on;
     plot(cal_60_end,sigma_eps_j(:,3),'-- k','LineWidth',3);
     plot(cal_60_end,sigma_eps_j(:,5),'-- k','LineWidth',3);
   hold off;
  title('\sigma_{\epsilon,{\iti}}');
  ax = gca;
  ax.FontSize = 20;
  ax.XLim = [1955 2020];


j = 2;
sigma_eps_j = sqrt(reshape(var_eps_unique_pct_17c(:,:,j),dnobs_60_end,6));

subplot(4,2,5);
   plot(cal_60_end,dp_disagg_60_end(:,j),'- k','LineWidth',1);
   tmp = char(labelvec_short_disagg(j));
  titstr =['(3) ' tmp]; 
  title(titstr);
  ax = gca;
  ax.FontSize = 20;
  ax.XLim = [1955 2020];
subplot(4,2,6);
  plot(cal_60_end,sigma_eps_j(:,1),'- k','LineWidth',3);
   hold on;
     plot(cal_60_end,sigma_eps_j(:,3),'-- k','LineWidth',3);
     plot(cal_60_end,sigma_eps_j(:,5),'-- k','LineWidth',3);
   hold off;
  title('\sigma_{\epsilon,{\iti}}');
  ax = gca;
  ax.FontSize = 20;
  ax.XLim = [1955 2020];
  
j = 7;
sigma_eps_j = sqrt(reshape(var_eps_unique_pct_17c(:,:,j),dnobs_60_end,6));

subplot(4,2,7);
   plot(cal_60_end,dp_disagg_60_end(:,j),'- k','LineWidth',1);
   tmp = char(labelvec_short_disagg(j));
  titstr =['(4) ' tmp]; 
  title(titstr);
  ax = gca;
  ax.FontSize = 20;
  ax.XLim = [1955 2020];
subplot(4,2,8);
  plot(cal_60_end,sigma_eps_j(:,1),'- k','LineWidth',3);
   hold on;
     plot(cal_60_end,sigma_eps_j(:,3),'-- k','LineWidth',3);
     plot(cal_60_end,sigma_eps_j(:,5),'-- k','LineWidth',3);
   hold off;
  title('\sigma_{\epsilon,{\iti}}');
  ax = gca;
  ax.FontSize = 20;
  ax.XLim = [1955 2020];

  