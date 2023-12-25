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
  matdir = 'matlab/mat/mucsv_tvp_17c_Q3/';
  
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
  ulabel = 'dp_agg_Q3';
  
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
  str_tmp = [matdir 'tau_mean_pct' ulabel]; load(str_tmp);
  
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
  
  
  figure;
  subplot(2,2,1);
   plot(cal_60_end,tau_mean_pct(:,1),'-- b','LineWidth',2);
   hold on;
     plot(cal_60_end,agg_tau_mean_pct(:,1),'- k','LineWidth',2);
     legend('UCSVO (PCE-all)','MUCSVO');
     legend('Location','NorthEast');
     hold off;
   titstr =['(a) Aggregate {\tau}_{\itt}' ]; 
   title(titstr);
   ax = gca;
   ax.FontSize = 20;
   ax.XLim = [1955 2020];
   
   subplot(2,2,2);
   plot(cal_60_end,sigma_dtau_common_mean_pct(:,4),'- k','LineWidth',3);
   hold on;
    plot(cal_60_end,sigma_dtau_common_mean_pct(:,3),'-- k','LineWidth',3);
    plot(cal_60_end,sigma_dtau_common_mean_pct(:,5),'-- k','LineWidth',3);
   hold off;
   titstr =['(b) {\sigma}_{{\Delta}{\tau} , {\itc} , {\itt}}' ]; 
   title(titstr);
   ax = gca;
   ax.FontSize = 20;
   ax.XLim = [1955 2020];
   
   subplot(2,2,3);
   plot(cal_60_end,sigma_eps_common_mean_pct(:,4),'- k','LineWidth',3);
   hold on;
    plot(cal_60_end,sigma_eps_common_mean_pct(:,3),'-- k','LineWidth',3);
    plot(cal_60_end,sigma_eps_common_mean_pct(:,5),'-- k','LineWidth',3);
   hold off;
   titstr =['(c) {\sigma}_{{\epsilon} , {\itc}  , {\itt}}' ]; 
   title(titstr);
   ax = gca;
   ax.FontSize = 20;
   ax.XLim = [1955 2020];
   
   subplot(2,2,4);
   plot(cal_60_end,scale_eps_common_mean_pct(:,4),'- k','LineWidth',3);
   hold on;
    plot(cal_60_end,scale_eps_common_mean_pct(:,3),'-- k','LineWidth',3);
    plot(cal_60_end,scale_eps_common_mean_pct(:,5),'-- k','LineWidth',3);
   hold off;
   titstr =['(d) {\its}_{{\itc} , {\itt}}' ]; 
   title(titstr);
   ylim([0 10]);
   ax = gca;
   ax.FontSize = 20;
   ax.XLim = [1955 2020];
    
  
  j = 15; % 'Fin. services & insurance'
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
   tmp = char(labelvec_disagg(j));
   titstr =['(a) ' tmp ]; 
   title(titstr);
   legend('Sectoral Inflation','Sectoral Trend');
   legend('Location','NorthEast');
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

  