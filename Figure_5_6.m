% Figures 5 and 6 
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
  
  % -- Some results for the 17-component model
  s = share_avg(ismpl_60_end==1,:);
  ifood = 5;
  ienergy = [7 10];
  icore = [1:4 6 8:9 11:17];
  s_energy = s(:,ienergy);
  s_energy = s_energy./repmat(sum(s_energy,2),1,2);
  s_core = s(:,icore);
  s_core = s_core./repmat(sum(s_core,2),1,14);
 
  % -- Mean Values of Parameters
  str_tmp = [matdir 'alpha_eps_mean' mlabel]; load(str_tmp);
  str_tmp = [matdir 'alpha_tau_mean' mlabel]; load(str_tmp);
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
  str_tmp = [matdir 'var_dtau_common_mean' mlabel]; load(str_tmp); 
  str_tmp = [matdir 'var_dtau_unique_mean' mlabel]; load(str_tmp);

  alpha_eps_17c = alpha_eps_mean;
  alpha_tau_17c = alpha_tau_mean;
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
  var_dtau_common_17c = var_dtau_common_mean;
  var_dtau_unique_17c = var_dtau_unique_mean;
  
  % ------------- Calculations for IRFs and weights used in Figure 5 -----------
  n_y = 17;
  alpha_eps = alpha_eps_17c;
  alpha_tau = alpha_tau_17c;
  sigma_eps_common = sqrt(var_eps_common_total_17c);
  sigma_dtau_common = sqrt(var_dtau_common_17c);
  sigma_eps_unique = sqrt(var_eps_unique_total_17c);
  sigma_dtau_unique = sqrt(var_dtau_unique_17c);
  
  agg_tau_decomp_total = zeros(dnobs_60_end,n_y,dnobs_60_end);
  for t = 1:dnobs_60_end;
      for i = 1:n_y;
          ydata = zeros(dnobs_60_end,n_y);
          ydata(t,i) = 1;
          [tmp1 tmp2 tmp3] = mfilt_agg(ydata,alpha_eps,alpha_tau,sigma_eps_common,sigma_dtau_common,sigma_eps_unique,sigma_dtau_unique,share_avg_60_end);
          tmp1 = tmp1(t:end);
          tmp2 = tmp2(t:end);
          tmp3 = tmp3(t:end);
          for j = 1:size(tmp1,1);
              agg_tau_decomp_total(t-1+j,i,j) = tmp1(j);
          end;
      end;
  end;
  
% Compute Sum of weights by series
wght = NaN*zeros(dnobs_60_end,n_y);
h = 4;
for i = 1:n_y;
    tmp = agg_tau_decomp_total(:,i,:);
    tmp = reshape(tmp,dnobs_60_end,dnobs_60_end);
    wght(:,i) = sum(tmp(:,1:h),2);
end;
tmp = sum(wght,2);
wght = wght./repmat(sum(wght,2),1,n_y);
 

  %------ Figure 5--------;
  j = 0;
  figure;
  for i = 1:n_y;
    j = j+1;
    subplot(6,3,j);
      titstr = labelvec_short_disagg(i);
      plot(cal_60_end,wght(:,i),'- b','LineWidth',3);
      hold on;
       plot(cal_60_end,share_avg_60_end(:,i),'-- r','LineWidth',3);
      hold off;
      title(titstr);
      ax = gca;
      ax.FontSize = 18;
      ylim([0 0.20]);
      ax.XLim = [1965 2015];
      if j == 1;
        legend('Weight','Share');
        legend('Location','NorthWest');
        %legend('Location','SouthEast');
      end;
      if j == 18;
          waitforbuttonpress;
          clf;
          j = 0;
          figure;
      end;
  end;
% error('tmp');
 
 
 %------ Figure 6--------;
 figure
  subplot(2,2,1);
      titstr = labelvec_short_disagg(i);
      plot(cal_60_end,sum(wght(:,icore),2),'- b','LineWidth',3);
      hold on;
       plot(cal_60_end,sum(share_avg_60_end(:,icore),2),'-- r','LineWidth',3);
      hold off;
      title('(a) Core','FontSize',18);
      ax = gca;
      ax.FontSize = 20;
      ylim([0 1.0]);
      ax.XLim = [1965 2015];
      legend('Weight','Share');
      legend('Location','SouthWest');
      
   subplot(2,2,2);
      titstr = labelvec_short_disagg(i);
      plot(cal_60_end,sum(wght(:,ifood),2),'- b','LineWidth',3);
      hold on;
       plot(cal_60_end,sum(share_avg_60_end(:,ifood),2),'-- r','LineWidth',3);
      hold off;
      title('(b) Food','FontSize',18);
      ax = gca;
      ax.FontSize = 20;
      ylim([0 1.0]);
      ax.XLim = [1965 2015];
    subplot(2,2,3);
      titstr = labelvec_short_disagg(i);
      plot(cal_60_end,sum(wght(:,ienergy),2),'- b','LineWidth',3);
      hold on;
       plot(cal_60_end,sum(share_avg_60_end(:,ienergy),2),'-- r','LineWidth',3);
      hold off;
      title('(c) Energy','FontSize',18);
      ax = gca;
      ax.FontSize = 20;
      ylim([0 1.0]);
      ax.XLim = [1965 2015];
  

  