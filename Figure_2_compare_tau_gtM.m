% comparing tau between model variants
% December 25th, 2023
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
  matdir_constvol_eps = 'matlab/mat/gtM-const_vol_e/';
  matdir_constvol_dtau = 'matlab/mat/gtM-const_vol_dtau/';

  % -- Read in Data --- 
  load_data = 1;  % 1 if reloading data from Excel, etc 
  mtoq_agg = 0;   % Temporal aggregation indicator of monthly to quarterly data
  nper = 12;
  pcomp_data_calendar_m_and_q_gt;
  
  % Data Series Used
  dp_agg = dp_agg_m;
  dp_agg_xfe = dp_agg_xfe_m;
  dp_agg_xe = dp_agg_xe_m;
  
  % dp_disagg = dp_disagg_q;
  % share_avg = share_avg_q;
  
  calvec = calvec_m;
  dnobs = dnobs_m;
  calds = calds_m;
  
  label_suffix = '_gtM';
  
  % Dates
  first_date = [2001 1];
  last_date = calds(end,:);
  ismpl_60_end = smpl(calvec,first_date,last_date,nper);
  cal_60_end = calvec(ismpl_60_end==1);
  dnobs_60_end = size(cal_60_end,1);
  
  dp_agg_60_end = dp_agg(ismpl_60_end==1,:);
  dp_agg_xe_60_end = dp_agg_xe(ismpl_60_end==1,:);
  dp_agg_xfe_60_end = dp_agg_xfe(ismpl_60_end==1,:);
  % dp_disagg_60_end = dp_disagg(ismpl_60_end==1,:);
  % share_avg_60_end = share_avg(ismpl_60_end==1,:);
  
  % -- Univariate Filtered and Smoothed Values
  ulabel = ['dp_agg' label_suffix];
  str_tmp = [matdir 'tau_mean_pct' ulabel]; load(str_tmp);
  tau_s_mean_pct_baseline = tau_mean_pct;
  str_tmp = [matdir 'var_dtau_mean_pct' ulabel]; load(str_tmp);
  var_dtau_mean_pct_baseline = var_dtau_mean_pct;
  str_tmp = [matdir 'var_eps_mean_pct' ulabel]; load(str_tmp);
  var_eps_mean_pct_baseline = var_eps_mean_pct;

  str_tmp = [matdir_nojumps 'tau_mean_pct' ulabel]; load(str_tmp);
  tau_s_mean_pct_nojumps = tau_mean_pct;
  str_tmp = [matdir_nojumps 'var_dtau_mean_pct' ulabel]; load(str_tmp);
  var_dtau_mean_pct_nojumps = var_dtau_mean_pct;
  str_tmp = [matdir_nojumps 'var_eps_mean_pct' ulabel]; load(str_tmp);
  var_eps_mean_pct_nojumps = var_eps_mean_pct;

  str_tmp = [matdir_constvol_eps 'tau_mean_pct' ulabel]; load(str_tmp);
  tau_s_mean_pct_constvol_eps = tau_mean_pct;
  str_tmp = [matdir_constvol_eps 'var_dtau_mean_pct' ulabel]; load(str_tmp);
  var_dtau_mean_pct_constvol_eps = var_dtau_mean_pct;
  str_tmp = [matdir_constvol_eps 'var_eps_mean_pct' ulabel]; load(str_tmp);
  var_eps_mean_pct_constvol_eps = var_eps_mean_pct;

  str_tmp = [matdir_constvol_dtau 'tau_mean_pct' ulabel]; load(str_tmp);
  tau_s_mean_pct_constvol_dtau = tau_mean_pct;
  str_tmp = [matdir_constvol_dtau 'var_dtau_mean_pct' ulabel]; load(str_tmp);
  var_dtau_mean_pct_constvol_dtau = var_dtau_mean_pct;
  str_tmp = [matdir_constvol_dtau 'var_eps_mean_pct' ulabel]; load(str_tmp);
  var_eps_mean_pct_constvol_dtau = var_eps_mean_pct;

  xdates_gt = [2000 2024];
  xticks_gt = 2000:2:2024; 
  font_size = 16; 
 
  %% --- Figure 2 comparing tau between model variants -----
  figure;
%   subplot(3,1,1);
  plot(cal_60_end,tau_s_mean_pct_baseline(:,1),'k-','LineWidth',4);
  hold on;
   plot(cal_60_end,tau_s_mean_pct_nojumps(:,1),'b-','LineWidth',2);
   plot(cal_60_end,tau_s_mean_pct_constvol_eps(:,1),'r:','LineWidth',2);
   plot(cal_60_end,tau_s_mean_pct_constvol_dtau(:,1),'m-','LineWidth',2);
  hold off;

  titstr = '(a) \tau_{{\it t}} ';
  title(titstr,'FontSize',18);
  legend('Baseline','No Outlier Adjustment','Constant Volatility \epsilon_{\itt}','Constant Volatility \Delta\tau_{\itt}');
  legend('Location','southoutside', 'Orientation','horizontal');
  xlim(xdates_gt);
  xticks(xticks_gt);
  ylim([-2 15]);
  ax = gca;
  ax.FontSize = font_size;
  
  
  %% --- Figure 3 comparing sigma_dtau between model variants -----
  figure;
%   subplot(3,1,2);
  plot(cal_60_end,var_dtau_mean_pct_baseline(:,1).^0.5,'k-','LineWidth',4);
  hold on;
   plot(cal_60_end,var_dtau_mean_pct_nojumps(:,1).^0.5,'b-','LineWidth',2);
   plot(cal_60_end,var_dtau_mean_pct_constvol_eps(:,1).^0.5,'r:','LineWidth',2);
   plot(cal_60_end,var_dtau_mean_pct_constvol_dtau(:,1).^0.5,'m-','LineWidth',2);
  hold off;

  titstr = '(b) \sigma_{\Delta\tau,{\itt}}';
  title(titstr,'FontSize',18);
  legend('Baseline','No Outlier Adjustment','Constant Volatility \epsilon_{\itt}','Constant Volatility \Delta\tau_{\itt}');
  legend('Location','southoutside', 'Orientation','horizontal');
  xlim(xdates_gt);
  xticks(xticks_gt);
  
  ax = gca;
%   ax.XTickLabel = []; 
  ax.FontSize = font_size;

  %% --- Figure 4 comparing sigma_eps between model variants -----
  figure;
%   subplot(3,1,3);
  plot(cal_60_end, var_eps_mean_pct_baseline(:,1).^0.5,'k-','LineWidth',4);
  hold on;
   plot(cal_60_end,var_eps_mean_pct_nojumps(:,1).^0.5,'b-','LineWidth',2);
   plot(cal_60_end,var_eps_mean_pct_constvol_eps(:,1).^0.5,'r:','LineWidth',2);
   plot(cal_60_end,var_eps_mean_pct_constvol_dtau(:,1).^0.5,'m-','LineWidth',2);
  hold off;

  titstr = '(c) \sigma_{\epsilon,{\itt}}';
  title(titstr,'FontSize',18);
  legend('Baseline','No Outlier Adjustment','Constant Volatility \epsilon_{\itt}','Constant Volatility \Delta\tau_{\itt}');
  legend('Location','southoutside', 'Orientation','horizontal');
  xlim(xdates_gt);
  xticks(xticks_gt);
  ax = gca;
  ax.FontSize = font_size;