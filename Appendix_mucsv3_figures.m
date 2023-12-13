% Summary Calculations 
% June 5, 2015
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
  %  Data Series Used
  dp_agg = dp_agg_q;
  dp_agg_xfe = dp_agg_xfe_q;
  dp_agg_xe = dp_agg_xe_q;
  dp_disagg = dp_disagg_3comp_q;
  share_avg = share_avg_3comp_q;
  share_avg_xfe = share_avg_3comp_xfe_q;
  share_avg_xe = share_avg_3comp_xe_q;
  calvec = calvec_q;
  dnobs = dnobs_q;
  calds = calds_q;
  nper = 4;
  mlabel = '_mucsv_tvp_3c_Q3';
  labelvec_disagg = labelvec_disagg_3comp;
  labelvec_short_disagg = labelvec_disagg_3comp;
  namevec_disagg = namevec_disagg_3comp;
  n_incl = n_incl_3comp;
  
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
  
  % Read UCSVO Estimates
  str_tmp = [matdir 'agg_tau_mean_pct' mlabel]; load(str_tmp);
  str_tmp = [matdir 'disagg_tau_mean' mlabel]; load(str_tmp);

  % -- Parameters
  str_tmp = [matdir 'alpha_eps_mean_pct' mlabel]; load(str_tmp);
  str_tmp = [matdir 'alpha_tau_mean_pct' mlabel]; load(str_tmp);
  str_tmp = [matdir 'sigma_eps_common_mean_pct' mlabel]; load(str_tmp);
  str_tmp = [matdir 'sigma_eps_unique_mean_pct' mlabel]; load(str_tmp);
  str_tmp = [matdir 'sigma_dtau_common_mean_pct' mlabel]; load(str_tmp); 
  str_tmp = [matdir 'sigma_dtau_unique_mean_pct' mlabel]; load(str_tmp);
  str_tmp = [matdir 'scale_eps_common_mean_pct' mlabel]; load(str_tmp);
  str_tmp = [matdir 'scale_eps_unique_mean_pct' mlabel]; load(str_tmp);
  

 for j = 1:3;
  alpha_eps_j = reshape(alpha_eps_mean_pct(:,:,j),dnobs_60_end,6);
  alpha_tau_j = reshape(alpha_tau_mean_pct(:,:,j),dnobs_60_end,6);
  sigma_eps_j = reshape(sigma_eps_unique_mean_pct(:,:,j),dnobs_60_end,6);
  sigma_dtau_j = reshape(sigma_dtau_unique_mean_pct(:,:,j),dnobs_60_end,6);
  scale_eps_j = reshape(scale_eps_unique_mean_pct(:,:,j),dnobs_60_end,6);
  
  figure;
  subplot(3,2,1);
   plot(cal_60_end,dp_disagg_60_end(:,j),'- k','LineWidth',1);
   hold on; 
     plot(cal_60_end,disagg_tau_mean(:,j),'-- r','LineWidth',3);
   hold off;
   tmp = char(labelvec_short_disagg(j));
   titstr =['(a) ' tmp ]; 
   title(titstr);
   legend('Sectoral Inflation','Sectoral Trend');
   legend('Location','SouthWest');
   ax = gca;
   ax.FontSize = 20;
   ax.XLim = [1955 2020];
   
   subplot(3,2,2);
   plot(cal_60_end,alpha_tau_j(:,4),'- k','LineWidth',3);
   hold on;
    plot(cal_60_end,alpha_tau_j(:,3),'-- k','LineWidth',3);
    plot(cal_60_end,alpha_tau_j(:,5),'-- k','LineWidth',3);
   hold off;
   tmp = char(labelvec_short_disagg(j));
   titstr =['(b) {\alpha}_{{\tau} , {\iti} , {\itt}}' ]; 
   title(titstr);
   ax = gca;
   ax.FontSize = 20;
   ax.XLim = [1955 2020];
   
   subplot(3,2,3);
   plot(cal_60_end,alpha_eps_j(:,4),'- k','LineWidth',3);
   hold on;
    plot(cal_60_end,alpha_eps_j(:,3),'-- k','LineWidth',3);
    plot(cal_60_end,alpha_eps_j(:,5),'-- k','LineWidth',3);
   hold off;
   tmp = char(labelvec_short_disagg(j));
   titstr =['(c) {\alpha}_{{\epsilon} , {\iti} , {\itt}}' ]; 
   title(titstr);
   ax = gca;
   ax.FontSize = 20;
   ax.XLim = [1955 2020];
   
   subplot(3,2,4);
   plot(cal_60_end,sigma_dtau_j(:,4),'- k','LineWidth',3);
   hold on;
    plot(cal_60_end,sigma_dtau_j(:,3),'-- k','LineWidth',3);
    plot(cal_60_end,sigma_dtau_j(:,5),'-- k','LineWidth',3);
   hold off;
   tmp = char(labelvec_short_disagg(j));
   titstr =['(d) {\sigma}_{{\Delta}{\tau} , {\iti} , {\itt}}' ]; 
   title(titstr);
   ax = gca;
   ax.FontSize = 20;
   ax.XLim = [1955 2020];
   
   subplot(3,2,5);
   plot(cal_60_end,sigma_eps_j(:,4),'- k','LineWidth',3);
   hold on;
    plot(cal_60_end,sigma_eps_j(:,3),'-- k','LineWidth',3);
    plot(cal_60_end,sigma_eps_j(:,5),'-- k','LineWidth',3);
   hold off;
   tmp = char(labelvec_short_disagg(j));
   titstr =['(e) {\sigma}_{{\epsilon} , {\iti} , {\itt}}' ]; 
   title(titstr);
   ax = gca;
   ax.FontSize = 20;
   ax.XLim = [1955 2020];
   
   subplot(3,2,6);
   plot(cal_60_end,scale_eps_j(:,4),'- k','LineWidth',3);
   hold on;
    plot(cal_60_end,scale_eps_j(:,3),'-- k','LineWidth',3);
    plot(cal_60_end,scale_eps_j(:,5),'-- k','LineWidth',3);
   hold off;
   tmp = char(labelvec_short_disagg(j));
   titstr =['(f) {\its}_{{\iti} , {\itt}}']; 
   title(titstr);
   ylim([0 10]);
   ax = gca;
   ax.FontSize = 20;
   ax.XLim = [1955 2020];

 end;