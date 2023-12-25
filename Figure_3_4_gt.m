% Summary Calculations 
% June 5, 2015
%


clear all;
small = 1.0e-10;
big = 1.0e+6;                  
rng(63761);
 
  % -- File Directories  
  outdir = 'matlab/out/';
  figdir = 'matlab/fig/';
  matdiruni = 'matlab/mat/gtM/';
  matdir = 'matlab/mat/gtM-12c/';
  
    % -- Read in Data --- 
  load_data = 1;  % 1 if reloading data from Excel, etc 
  mtoq_agg = 0;   % Temporal aggregation indicator of monthly to quarterly data
  pcomp_data_calendar_m_and_q_gt;
  nper = 12;
  
  % Data Series Used
  dp_agg = dp_agg_m;
%   dp_agg_xfe = dp_agg_xfe_q;
%   dp_agg_xe = dp_agg_xe_q;
  dp_disagg = dp_disagg_m;
  share_avg = share_avg_m;
  
  % Univariate calendars
  calvec_uni = calvec_m; 
  calds_uni = calds_m;
  
  % MUSCVO calendars
  calvec = calvec_disagg_m;
  dnobs = dnobs_disagg_m;
  calds = calds_disagg_m;
  
  mlabel = '_mucsv_tvp_12c_gt_M';
  ulabel = 'dp_agg_gtM';
  
  % Dates
  first_date = [2011 1];
  last_date = calds(end,:);
  ismpl_60_end = smpl(calvec,first_date,last_date, nper);
  cal_60_end = calvec(ismpl_60_end==1);
  dnobs_60_end = size(cal_60_end,1);
  
  dp_agg_60_end = dp_agg(ismpl_60_end==1,:);
%   dp_agg_xe_60_end = dp_agg_xe(ismpl_60_end==1,:);
%   dp_agg_xfe_60_end = dp_agg_xfe(ismpl_60_end==1,:);
  dp_disagg_60_end = dp_disagg(ismpl_60_end==1,:);
  share_avg_60_end = share_avg(ismpl_60_end==1,:);
  
  % Read UCSVO estimates
  str_tmp = [matdiruni 'tau_mean_pct' ulabel]; load(str_tmp);
  
  % Read MUSCVO estimates
  str_tmp = [matdir 'agg_tau_mean_pct' mlabel]; load(str_tmp);
  str_tmp = [matdir 'disagg_tau_mean' mlabel]; load(str_tmp);

  % -- Parameters of MUCSVO
  str_tmp = [matdir 'alpha_eps_mean_pct' mlabel]; load(str_tmp);
  str_tmp = [matdir 'alpha_tau_mean_pct' mlabel]; load(str_tmp);
  str_tmp = [matdir 'sigma_eps_common_mean_pct' mlabel]; load(str_tmp);
  str_tmp = [matdir 'sigma_eps_unique_mean_pct' mlabel]; load(str_tmp);
  str_tmp = [matdir 'sigma_dtau_common_mean_pct' mlabel]; load(str_tmp); 
  str_tmp = [matdir 'sigma_dtau_unique_mean_pct' mlabel]; load(str_tmp);
  str_tmp = [matdir 'scale_eps_common_mean_pct' mlabel]; load(str_tmp);
  str_tmp = [matdir 'scale_eps_unique_mean_pct' mlabel]; load(str_tmp);
  
  %% Aggregate and common unobserved component
  figure;
  subplot(2,2,1);
   plot(calvec_uni(2:end),tau_mean_pct(:,1),'- b','LineWidth',1.5);
   hold on;
     plot(cal_60_end,agg_tau_mean_pct(:,1),'- k','LineWidth',3);
     legend('UCSVO (CPI-all)','MUCSVO');
     legend('Location','NorthEast');
     hold off;
   titstr =['(a) Aggregate {\tau}_{\itt}' ]; 
   title(titstr);
   ax = gca;
   ax.FontSize = 16;
   ax.XLim = [2000 2024];
   
   % Volatility of dtau 
   subplot(2,2,2);
   plot(cal_60_end,sigma_dtau_common_mean_pct(:,4),'- k','LineWidth',3);
   hold on;
    plot(cal_60_end,sigma_dtau_common_mean_pct(:,3),'-- k','LineWidth',2);
    plot(cal_60_end,sigma_dtau_common_mean_pct(:,5),'-- k','LineWidth',2);
   hold off;
   titstr =['(b) {\sigma}_{{\Delta}{\tau} , {\itc} , {\itt}}' ]; 
   title(titstr);
   ax = gca;
   ax.FontSize = 16;
   ax.XLim = [2010 2024];
   
   % Volatility of unobserved common component
   subplot(2,2,3);
   plot(cal_60_end,sigma_eps_common_mean_pct(:,4),'- k','LineWidth',3);
   hold on;
    plot(cal_60_end,sigma_eps_common_mean_pct(:,3),'-- k','LineWidth',2);
    plot(cal_60_end,sigma_eps_common_mean_pct(:,5),'-- k','LineWidth',2);
   hold off;
   titstr =['(c) {\sigma}_{{\epsilon} , {\itc}  , {\itt}}' ]; 
   title(titstr);
   ax = gca;
   ax.FontSize = 16;
   ax.XLim = [2010 2024];
   
   % s_t of unobserved common component
   subplot(2,2,4);
   plot(cal_60_end,scale_eps_common_mean_pct(:,4),'- k','LineWidth',3);
   hold on;
    plot(cal_60_end,scale_eps_common_mean_pct(:,3),'-- k','LineWidth',2);
    plot(cal_60_end,scale_eps_common_mean_pct(:,5),'-- k','LineWidth',2);
   hold off;
   titstr =['(d) {\its}_{{\itc} , {\itt}}' ]; 
   title(titstr);
   ylim([0 10]);
   ax = gca;
   ax.FontSize = 16;
   ax.XLim = [2010 2024];
    
  
  %%  --- Specific sector plot for Food and Non-alcoholic beverages ---
  
  j = 1; % Sector j
  alpha_eps_j = reshape(alpha_eps_mean_pct(:,:,j),dnobs_60_end,6);
  alpha_tau_j = reshape(alpha_tau_mean_pct(:,:,j),dnobs_60_end,6);
  sigma_eps_j = reshape(sigma_eps_unique_mean_pct(:,:,j),dnobs_60_end,6);
  sigma_dtau_j = reshape(sigma_dtau_unique_mean_pct(:,:,j),dnobs_60_end,6);
  scale_eps_j = reshape(scale_eps_unique_mean_pct(:,:,j),dnobs_60_end,6);
  
  figure;
  subplot(3,2,1);
   plot(cal_60_end,dp_disagg_60_end(:,j),'- k','LineWidth',1);
   hold on; 
     plot(cal_60_end,disagg_tau_mean(:,j),'-- r','LineWidth',2);
   hold off;
   % tmp = char(labelvec_disagg(j));
   tmp = 'Food and Non-Alcoholic Beverages';
   titstr =['(a) ' tmp ]; 
   title(titstr);
   legend('Sectoral Inflation','Sectoral Trend');
   legend('Location','NorthEast');
   ax = gca;
   ax.FontSize = 16;
   ax.XLim = [2010 2024];
   
   subplot(3,2,2);
   plot(cal_60_end,alpha_tau_j(:,4),'- k','LineWidth',3);
   hold on;
    plot(cal_60_end,alpha_tau_j(:,3),'-- k','LineWidth',2);
    plot(cal_60_end,alpha_tau_j(:,5),'-- k','LineWidth',2);
   hold off;
%    tmp = char(labelvec_short_disagg(j));
   titstr =['(b) {\alpha}_{{\tau} , {\iti} , {\itt}}' ]; 
   title(titstr);
   ax = gca;
   ax.FontSize = 16;
   ax.XLim = [2010 2024];
   ax.YLim = [0.1 1.2];
   
   subplot(3,2,3);
   plot(cal_60_end,alpha_eps_j(:,4),'- k','LineWidth',3);
   hold on;
    plot(cal_60_end,alpha_eps_j(:,3),'-- k','LineWidth',2);
    plot(cal_60_end,alpha_eps_j(:,5),'-- k','LineWidth',2);
   hold off;
%    tmp = char(labelvec_short_disagg(j));
   titstr =['(c) {\alpha}_{{\epsilon} , {\iti} , {\itt}}' ]; 
   title(titstr);
   ax = gca;
   ax.FontSize = 16;
   ax.XLim = [2010 2024];
   ax.YLim = [-0.6 0.6];
   
   subplot(3,2,4);
   plot(cal_60_end,sigma_dtau_j(:,4),'- k','LineWidth',3);
   hold on;
    plot(cal_60_end,sigma_dtau_j(:,3),'-- k','LineWidth',2);
    plot(cal_60_end,sigma_dtau_j(:,5),'-- k','LineWidth',2);
   hold off;
%    tmp = char(labelvec_short_disagg(j));
   titstr =['(d) {\sigma}_{{\Delta}{\tau} , {\iti} , {\itt}}' ]; 
   title(titstr);
   ax = gca;
   ax.FontSize = 16;
   ax.XLim = [2010 2024];
   ax.YLim = [0 1];
   
   subplot(3,2,5);
   plot(cal_60_end,sigma_eps_j(:,4),'- k','LineWidth',3);
   hold on;
    plot(cal_60_end,sigma_eps_j(:,3),'-- k','LineWidth',2);
    plot(cal_60_end,sigma_eps_j(:,5),'-- k','LineWidth',2);
   hold off;
%    tmp = char(labelvec_short_disagg(j));
   titstr =['(e) {\sigma}_{{\epsilon} , {\iti} , {\itt}}' ]; 
   title(titstr);
   ax = gca;
   ax.FontSize = 16;
   ax.XLim = [2010 2024];
   ax.YLim = [4 16];
   
   subplot(3,2,6);
   plot(cal_60_end,scale_eps_j(:,4),'- k','LineWidth',3);
   hold on;
    plot(cal_60_end,scale_eps_j(:,3),'-- k','LineWidth',2);
    plot(cal_60_end,scale_eps_j(:,5),'-- k','LineWidth',2);
   hold off;
%    tmp = char(labelvec_short_disagg(j));
   titstr =['(f) {\its}_{{\iti} , {\itt}}']; 
   title(titstr);
   ylim([0 10]);
   ax = gca;
   ax.FontSize = 16;
   ax.XLim = [2010 2024];

  