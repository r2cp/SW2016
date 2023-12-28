% Figures 1 and 2
% November 4
%

% clear all;
small = 1.0e-10;
big = 1.0e+6;                  
rng(63761);
 
 % -- File Directories  
  outdir = 'matlab/out/';
  figdir = 'matlab/fig/';
  matdir = 'matlab/mat/gtM-const_vol_dtau/';
  
  % -- Read in Data --- 
  load_data = 1;  % 1 if reloading data from Excel, etc 
  mtoq_agg = 0;   % Temporal aggregation indicator of monthly to quarterly data
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
  nper = 12;
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
  tau_s_mean_pct_h = tau_mean_pct;
  str_tmp = [matdir 'var_dtau_mean_pct' ulabel]; load(str_tmp);
  str_tmp = [matdir 'vartotal_eps_mean_pct' ulabel]; load(str_tmp);
  str_tmp = [matdir 'var_eps_mean_pct' ulabel]; load(str_tmp);
  str_tmp = [matdir 'scl_eps_pct' ulabel]; load(str_tmp);
  var_eps_mean_pct_h=var_eps_mean_pct;
  vartotal_eps_mean_pct_h=vartotal_eps_mean_pct;
  var_dtau_mean_pct_h = var_dtau_mean_pct;
  scl_eps_pct_h = scl_eps_pct;
   
  ulabel = ['dp_xe' label_suffix];
  str_tmp = [matdir 'tau_mean_pct' ulabel]; load(str_tmp);
  tau_s_mean_pct_xe = tau_mean_pct;
  str_tmp = [matdir 'var_dtau_mean_pct' ulabel]; load(str_tmp);
  str_tmp = [matdir 'vartotal_eps_mean_pct' ulabel]; load(str_tmp);
  str_tmp = [matdir 'var_eps_mean_pct' ulabel]; load(str_tmp);
  str_tmp = [matdir 'scl_eps_pct' ulabel]; load(str_tmp);
  var_eps_mean_pct_xe=var_eps_mean_pct;
  vartotal_eps_mean_pct_xe=vartotal_eps_mean_pct;
  var_dtau_mean_pct_xe = var_dtau_mean_pct;
  scl_eps_pct_xe = scl_eps_pct;
  
  ulabel = ['dp_xfe' label_suffix];
  str_tmp = [matdir 'tau_mean_pct' ulabel]; load(str_tmp);
  tau_s_mean_pct_xfe = tau_mean_pct;
  str_tmp = [matdir 'var_dtau_mean_pct' ulabel]; load(str_tmp);
  str_tmp = [matdir 'vartotal_eps_mean_pct' ulabel]; load(str_tmp);
  str_tmp = [matdir 'var_eps_mean_pct' ulabel]; load(str_tmp);
  str_tmp = [matdir 'scl_eps_pct' ulabel]; load(str_tmp);
  var_eps_mean_pct_xfe=var_eps_mean_pct;
  vartotal_eps_mean_pct_xfe=vartotal_eps_mean_pct;
  var_dtau_mean_pct_xfe = var_dtau_mean_pct;
  scl_eps_pct_xfe = scl_eps_pct;
  
  xdates_gt = [2000 2024];

  % -- Figure 1: 4 Panel Graph
  figure;
  subplot(2,2,1);
  plot(cal_60_end,dp_agg_60_end,'- k','LineWidth',2);
  titstr = '(a) Headline CPI';
  title(titstr,'FontSize',18);
  xlim(xdates_gt);
  ylim([-10 15]);
  ax = gca;
  ax.FontSize = 20;
  
  subplot(2,2,2);
  plot(cal_60_end,dp_agg_xe_60_end,'- k','LineWidth',2);
  titstr = '(b) CPIxE Inflation ';
  title(titstr,'FontSize',18);
  xlim(xdates_gt);
  ylim([-10 15]);
  ax = gca;
  ax.FontSize = 20;
  
  subplot(2,2,3);
  plot(cal_60_end,dp_agg_xfe_60_end,'- k','LineWidth',2);
  titstr = '(c) CPIxFE Inflation ';
  title(titstr,'FontSize',18);
  xlim(xdates_gt);
  ylim([-10 15]);
  ax = gca;
  ax.FontSize = 20;
  
  
  
%   % MADeviation reported in paper 
%  [cal_60_end/1000 tau_s_h tau_s_xe tau_s_xfe]
%  tmp = abs(tau_s_h-tau_s_xfe);
%  mean(tmp)
%  [aa,ii]=max(tmp);
%  cal_60_end(ii)
%  aa
%  tmp = abs(tau_s_xe-tau_s_xfe);
%  mean(tmp)
%  [aa,ii]=max(tmp);
%  cal_60_end(ii)
%  aa

  
 % -- Figure 2 -----
figure;

xticks_gt = 2000:2:2023; 
font_size = 16; 

subplot(2,2,1); 
  plot(cal_60_end,tau_s_mean_pct_h(:,1),'- b','LineWidth',2);
  hold on;
   plot(cal_60_end,tau_s_mean_pct_xe(:,1),'-- k','LineWidth',2);
   plot(cal_60_end,tau_s_mean_pct_xfe(:,1),': r','LineWidth',2);
  hold off;
  titstr = '(a) \tau_{{\itt}} ';
  title(titstr,'FontSize',18);
  legend('CPI','CPIxE','CPIxFE');
  % legend('CPI','CPIxE');
  legend('Location','southoutside', 'Orientation','horizontal');
  xlim(xdates_gt);
  ylim([-2 15]);
  xticks(xticks_gt);
  ax = gca;
  ax.FontSize = font_size;
 
  subplot(2,2,2);
  plot(cal_60_end,sqrt(var_dtau_mean_pct_h(:,4)),'- b','LineWidth',2);
  hold on;
    plot(cal_60_end,sqrt(var_dtau_mean_pct_xe(:,4)),'-- k','LineWidth',2);
    plot(cal_60_end,sqrt(var_dtau_mean_pct_xfe(:,4)),': r','LineWidth',2);  
  hold off;
  titstr = '(b) \sigma_{\Delta\tau,{\itt}}';
  title(titstr,'FontSize',18);
  xlim(xdates_gt); 
  xticks(xticks_gt);
  set(gca,'fontsize',font_size)

  
subplot(2,2,3);
  plot(cal_60_end,sqrt(var_eps_mean_pct_h(:,4)),'- b','LineWidth',2);
  hold on;
    plot(cal_60_end,sqrt(var_eps_mean_pct_xe(:,4)),'-- k','LineWidth',2);
    plot(cal_60_end,sqrt(var_eps_mean_pct_xfe(:,4)),': r','LineWidth',2);  
  hold off;
  titstr = '(c) \sigma_{\epsilon,{\itt}}';
  title(titstr,'FontSize',18);
  xlim(xdates_gt);
  xticks(xticks_gt);
  ylim([0 6]);
  set(gca,'fontsize',font_size)
  
 subplot(2,2,4);
  plot(cal_60_end,scl_eps_pct_h(:,4),'- b','LineWidth',2);
  hold on;
    plot(cal_60_end,scl_eps_pct_xe(:,4),'-- k','LineWidth',2);
    plot(cal_60_end,scl_eps_pct_xfe(:,4),': r','LineWidth',2);  
  hold off;
  titstr = '(d) {\its}_{{\itt}}';
  title(titstr,'FontSize',18);
  xlim(xdates_gt);
  xticks(xticks_gt);
  set(gca,'fontsize',font_size)

 
 