% Figure 8
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
  label_suffix = '_Q3';
   
  % Dates
  first_date = [1960 1];
  last_date = calds(end,:);
  ismpl_60_end = smpl(calvec,first_date,last_date,nper);
  cal_60_end = calvec(ismpl_60_end==1);
  dnobs_60_end = size(cal_60_end,1);
  
  dp_agg_60_end = dp_agg(ismpl_60_end==1,:);
  dp_agg_xe_60_end = dp_agg_xe(ismpl_60_end==1,:);
  dp_agg_xfe_60_end = dp_agg_xfe(ismpl_60_end==1,:);
  dp_disagg_60_end = dp_disagg(ismpl_60_end==1,:);
  share_avg_60_end = share_avg(ismpl_60_end==1,:);
  
  % -- Univariate Filtered and Smoothed Values
  ulabel = ['dp_agg' label_suffix];
  str_tmp = [matdir 'tau_mean_pct' ulabel]; load(str_tmp);
  tau_s_h = tau_mean_pct(:,1);
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
  tau_s_xe = tau_mean_pct(:,1);
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
  tau_s_xfe = tau_mean_pct(:,1);
  tau_s_mean_pct_xfe = tau_mean_pct;
  str_tmp = [matdir 'var_dtau_mean_pct' ulabel]; load(str_tmp);
  str_tmp = [matdir 'vartotal_eps_mean_pct' ulabel]; load(str_tmp);
  str_tmp = [matdir 'var_eps_mean_pct' ulabel]; load(str_tmp);
  str_tmp = [matdir 'scl_eps_pct' ulabel]; load(str_tmp);
  var_eps_mean_pct_xfe=var_eps_mean_pct;
  vartotal_eps_mean_pct_xfe=vartotal_eps_mean_pct;
  var_dtau_mean_pct_xfe = var_dtau_mean_pct;
  scl_eps_pct_xfe = scl_eps_pct;
  
  % -- Some results for the 17-component model
  mlabel = ['_mucsv_tvp_17c' label_suffix];
  str_tmp = [matdir 'agg_tau_mean_pct' mlabel]; load(str_tmp);
  tau_s_17c = agg_tau_mean_pct(:,1);
  tau_s_mean_pct_17c = agg_tau_mean_pct;
  
  mlabel = ['_mucsv_tvp_17c' label_suffix];
  str_tmp = [matdir 'agg_tau_xfe_mean_pct' mlabel]; load(str_tmp);
  tau_xfe_s_17c = agg_tau_xfe_mean_pct(:,1);
  tau_xfe_s_mean_pct_17c = agg_tau_xfe_mean_pct;
  
  mlabel = ['_mucsv_tvp_17c' label_suffix];
  str_tmp = [matdir 'agg_tau_xe_mean_pct' mlabel]; load(str_tmp);
  tau_xe_s_17c = agg_tau_xe_mean_pct(:,1);
  tau_xe_s_mean_pct_17c = agg_tau_xe_mean_pct;
 
  % -- Some results for the 3-component model
  mlabel = ['_mucsv_tvp_3c' label_suffix];
  str_tmp = [matdir 'agg_tau_mean_pct' mlabel]; load(str_tmp);
  tau_s_3c = agg_tau_mean_pct(:,1);
  tau_s_mean_pct_3c = agg_tau_mean_pct;
  
  mlabel = ['_mucsv_tvp_3c' label_suffix];
  str_tmp = [matdir 'agg_tau_xfe_mean_pct' mlabel]; load(str_tmp);
  tau_xfe_s_3c = agg_tau_xfe_mean_pct(:,1);
  tau_xfe_s_mean_pct_3c = agg_tau_xfe_mean_pct;
  
  mlabel = ['_mucsv_tvp_3c' label_suffix];
  str_tmp = [matdir 'agg_tau_xe_mean_pct' mlabel]; load(str_tmp);
  tau_xe_s_3c = agg_tau_xe_mean_pct(:,1);
  tau_xe_s_mean_pct_3c = agg_tau_xe_mean_pct;
  
  
  i1 = 2;
  i2 = 6;
  % -- Figure 8: 4 Panel Graph
  figure;
  subplot(2,2,1);
  plot(cal_60_end,tau_s_mean_pct_h(:,i1),'- k','LineWidth',1);
  hold on;
   plot(cal_60_end,tau_s_mean_pct_17c(:,i1),'-- b','LineWidth',2);
   legend('Univariate','Multivariate (17 Components)');
   legend('Location','NorthEast');
   plot(cal_60_end,tau_s_mean_pct_h(:,i2),'- k','LineWidth',1);
   plot(cal_60_end,tau_s_mean_pct_17c(:,i2),'-- b','LineWidth',2);
  hold off
  titstr = '(a) {\tau}_{\itt}^{PCE-All}';
  title(titstr,'FontSize',18);
  xlim([1955 2020]);
  ylim([-2 14]);
  ax = gca;
  ax.FontSize = 20;
  
 
  subplot(2,2,2)
  plot(cal_60_end,tau_s_mean_pct_xe(:,2),'- k','LineWidth',1);
  hold on;
   plot(cal_60_end,tau_xe_s_mean_pct_17c(:,2),'-- b','LineWidth',2);
   legend('Univariate','Multivariate (17 Components)');
   legend('Location','NorthEast');
   plot(cal_60_end,tau_s_mean_pct_xe(:,6),'- k','LineWidth',1);
   plot(cal_60_end,tau_xe_s_mean_pct_17c(:,6),'-- b','LineWidth',2);
  hold off
  titstr = '(b) {\tau}_{\itt}^{PCExE}';
  title(titstr,'FontSize',18);
  xlim([1955 2020]);
  ylim([-2 14]);
  ax = gca;
  ax.FontSize = 20;
 
  subplot(2,2,3);
  plot(cal_60_end,tau_s_mean_pct_xfe(:,2),'- k','LineWidth',1);
  hold on;
   plot(cal_60_end,tau_xfe_s_mean_pct_17c(:,2),'-- b','LineWidth',2);
   legend('Univariate','Multivariate (17 Components)');
   legend('Location','NorthEast');
   plot(cal_60_end,tau_s_mean_pct_xfe(:,6),'- k','LineWidth',1);
   plot(cal_60_end,tau_xfe_s_mean_pct_17c(:,6),'-- b','LineWidth',2);
  hold off
  titstr = '(c) {\tau}_{\itt}^{PCExFE}';
  title(titstr,'FontSize',18);
  xlim([1955 2020]);
  ylim([-2 14]);
  ax = gca;
  ax.FontSize = 20;
  
  subplot(2,2,4);
  plot(cal_60_end,tau_s_mean_pct_3c(:,i1),'- k','LineWidth',1);
  hold on;
   plot(cal_60_end,tau_s_mean_pct_17c(:,i1),'-- b','LineWidth',2);
   legend('Multivariate (3 Components)','Multivariate (17 Components)');
   legend('Location','NorthEast');
   plot(cal_60_end,tau_s_mean_pct_3c(:,i2),'- k','LineWidth',1);
   plot(cal_60_end,tau_s_mean_pct_17c(:,i2),'-- b','LineWidth',2);
  hold off
  titstr = '(d) {\tau}_{\itt}^{PCE-All}';
  title(titstr,'FontSize',18);
  xlim([1955 2020]);
  ylim([-2 14]);
  ax = gca;
  ax.FontSize = 20;
  
  
  % Results in panel A of Table 2
  % Compute Average width of error bands
  ismpl1 = smpl(cal_60_end,[1960 1],[1989 4],4);
  ismpl2 = smpl(cal_60_end,[1990 1],last_date,4);
  i1 = 2; i2 = 6;
  uni_h_90 = tau_s_mean_pct_h(:,i2) - tau_s_mean_pct_h(:,i1);
  uni_xe_90 = tau_s_mean_pct_xe(:,i2) - tau_s_mean_pct_xe(:,i1);
  uni_xfe_90 = tau_s_mean_pct_xfe(:,i2) - tau_s_mean_pct_xfe(:,i1);
  m17_h_90 = tau_s_mean_pct_17c(:,i2)-tau_s_mean_pct_17c(:,i1);
  m17_xe_90 = tau_xe_s_mean_pct_17c(:,i2)-tau_xe_s_mean_pct_17c(:,i1);
  m17_xfe_90 = tau_xfe_s_mean_pct_17c(:,i2)-tau_xfe_s_mean_pct_17c(:,i1);
  m3_h_90 = tau_s_mean_pct_3c(:,i2)-tau_s_mean_pct_3c(:,i1);
  m3_xe_90 = tau_xe_s_mean_pct_3c(:,i2)-tau_xe_s_mean_pct_3c(:,i1);
  m3_xfe_90 = tau_xfe_s_mean_pct_3c(:,i2)-tau_xfe_s_mean_pct_3c(:,i1);
  i1 = 3; i2 = 5;
  uni_h_67 = tau_s_mean_pct_h(:,i2) - tau_s_mean_pct_h(:,i1);
  uni_xe_67 = tau_s_mean_pct_xe(:,i2) - tau_s_mean_pct_xe(:,i1);
  uni_xfe_67 = tau_s_mean_pct_xfe(:,i2) - tau_s_mean_pct_xfe(:,i1);
  m17_h_67 = tau_s_mean_pct_17c(:,i2)-tau_s_mean_pct_17c(:,i1);
  m17_xe_67 = tau_xe_s_mean_pct_17c(:,i2)-tau_xe_s_mean_pct_17c(:,i1);
  m17_xfe_67 = tau_xfe_s_mean_pct_17c(:,i2)-tau_xfe_s_mean_pct_17c(:,i1);
  m3_h_67 = tau_s_mean_pct_3c(:,i2)-tau_s_mean_pct_3c(:,i1);
  m3_xe_67 = tau_xe_s_mean_pct_3c(:,i2)-tau_xe_s_mean_pct_3c(:,i1);
  m3_xfe_67 = tau_xfe_s_mean_pct_3c(:,i2)-tau_xfe_s_mean_pct_3c(:,i1);
  
  % Output 
  outfile_name = [outdir 'Table_ErrorBand_Width' label_suffix '.out'];
  fileID = fopen(outfile_name,'w');
  fprintf(fileID,'Width of 67 percent bands \n');
  fprintf(fileID,'%4.2f,',mean(uni_h_67(ismpl1==1)));
  fprintf(fileID,'%4.2f,',mean(uni_h_67(ismpl2==1)));
  fprintf(fileID,'%4.2f,',mean(m3_h_67(ismpl1==1)));
  fprintf(fileID,'%4.2f,',mean(m3_h_67(ismpl2==1)));
  fprintf(fileID,'%4.2f,',mean(m17_h_67(ismpl1==1)));
  fprintf(fileID,'%4.2f \n',mean(m17_h_67(ismpl2==1)));
  fprintf(fileID,'%4.2f,',mean(uni_xe_67(ismpl1==1)));
  fprintf(fileID,'%4.2f,',mean(uni_xe_67(ismpl2==1)));
  fprintf(fileID,'%4.2f,',mean(m3_xe_67(ismpl1==1)));
  fprintf(fileID,'%4.2f,',mean(m3_xe_67(ismpl2==1)));
  fprintf(fileID,'%4.2f,',mean(m17_xe_67(ismpl1==1)));
  fprintf(fileID,'%4.2f \n',mean(m17_xe_67(ismpl2==1)));
  fprintf(fileID,'%4.2f,',mean(uni_xfe_67(ismpl1==1)));
  fprintf(fileID,'%4.2f,',mean(uni_xfe_67(ismpl2==1)));
  fprintf(fileID,'%4.2f,',mean(m3_xfe_67(ismpl1==1)));
  fprintf(fileID,'%4.2f,',mean(m3_xfe_67(ismpl2==1)));
  fprintf(fileID,'%4.2f,',mean(m17_xfe_67(ismpl1==1)));
  fprintf(fileID,'%4.2f \n',mean(m17_xfe_67(ismpl2==1)));
  
  fprintf(fileID,'Width of 90 percent bands \n');
  fprintf(fileID,'%4.2f,',mean(uni_h_90(ismpl1==1)));
  fprintf(fileID,'%4.2f,',mean(uni_h_90(ismpl2==1)));
  fprintf(fileID,'%4.2f,',mean(m3_h_90(ismpl1==1)));
  fprintf(fileID,'%4.2f,',mean(m3_h_90(ismpl2==1)));
  fprintf(fileID,'%4.2f,',mean(m17_h_90(ismpl1==1)));
  fprintf(fileID,'%4.2f \n',mean(m17_h_90(ismpl2==1)));
  fprintf(fileID,'%4.2f,',mean(uni_xe_90(ismpl1==1)));
  fprintf(fileID,'%4.2f,',mean(uni_xe_90(ismpl2==1)));
  fprintf(fileID,'%4.2f,',mean(m3_xe_90(ismpl1==1)));
  fprintf(fileID,'%4.2f,',mean(m3_xe_90(ismpl2==1)));
  fprintf(fileID,'%4.2f,',mean(m17_xe_90(ismpl1==1)));
  fprintf(fileID,'%4.2f \n',mean(m17_xe_90(ismpl2==1)));
  fprintf(fileID,'%4.2f,',mean(uni_xfe_90(ismpl1==1)));
  fprintf(fileID,'%4.2f,',mean(uni_xfe_90(ismpl2==1)));
  fprintf(fileID,'%4.2f,',mean(m3_xfe_90(ismpl1==1)));
  fprintf(fileID,'%4.2f,',mean(m3_xfe_90(ismpl2==1)));
  fprintf(fileID,'%4.2f,',mean(m17_xfe_90(ismpl1==1)));
  fprintf(fileID,'%4.2f \n',mean(m17_xfe_90(ismpl2==1)));
