% ucsv_ep_poos.m 
% 9/28/2015, MWW
% poos analysis
%


clear all;
small = 1.0e-10;
big = 1.0e+6;                  
rng(663436);
 
  % -- File Directories  
  outdir = '/Users/mwatson/Dropbox/p_comp/revision/matlab/out/';
  figdir = '/Users/mwatson/Dropbox/p_comp/revision/matlab/fig/';
  matdir = 'matlab/mat/gtM-const_vol_e-poos/';

  
  % -- Read in Data --- 
  load_data = 1;  % 1 if reloading data from Excel, etc 
  mtoq_agg = 0;   % Temporal aggregation indicator of monthly to quarterly data
  pcomp_data_calendar_m_and_q_gt;  
  
  % Data Series Used
  dp_agg = dp_agg_m;
  dp_agg_xfe = dp_agg_xfe_m;
  dp_agg_xe = dp_agg_xe_m;
  dp_disagg = dp_disagg_m;
  calvec = calvec_m;
  dnobs = dnobs_m;
  calds = calds_m;
  nper = 12;
  outlabel = 'gtM';
   
  labvec = {'Headline Inflation';'Core XFE';'Core XE'};
  namevec = {'dp_agg';'dp_xfe';'dp_xe'};
  dp = [dp_agg dp_agg_xfe dp_agg_xe];
%   dp = [dp_agg];

  % Dates
  first_date = [2001 1];
  last_date = calds(end,:);
  ismpl = smpl(calvec,first_date,last_date,nper);
  calvec_ismpl = calvec(ismpl==1);
  calds_ismpl = calds(ismpl==1,:);
  dnobs_ismpl = size(calvec_ismpl,1);
  str_tmp = [matdir 'calvec_ismpl_ucsv_poos']; save(str_tmp,'calvec_ismpl');
  
  % First POOS Date
  first_date_poos = [2015 1];
  tmp = calds_ismpl == repmat(first_date_poos,dnobs_ismpl,1);
  [a,i_first_poos] = max(sum(tmp,2));  % Index in calvec_ismpl of first poos period
  
  % Parameters for UCSV Draws
  % Draws for pre-poos period
  n_burnin = 1000;     % Discarded Draws
  n_draws_save = 1000;  % Number of Draws to Save
  k_draws = 5;          % Save results every k_draws
  
  % Parameters for scale mixture of epsilon component
  scl_eps_vec = [1; linspace(2.0,10.0,9)'];
  n_scl_eps = length(scl_eps_vec);
  ps_mean = 1 - 1/(2*nper);              % Outlier every 4 years
  ps_prior_obs = nper*10;                % Sample size of 10 years for prior
  ps_prior_a = ps_mean*ps_prior_obs;     % "alpha" in beta prior
  ps_prior_b = (1-ps_mean)*ps_prior_obs; % "beta" in beta prior
  ps_prior = [ps_prior_a ps_prior_b];
  
  % -- Parameters for RW Innovation Variance -- 
  % Gvalues .. for evolution of standard deviations over a year (nper periods)
  ng = 5;      % Number of grid points for approximate uniform prior
  g_dtau_min = 0.001;
  g_dtau_max = 0.20;
  g_dtau_values = linspace(g_dtau_min,g_dtau_max,ng)';
  p_g_dtau_values = ones(ng,1)/ng;
  g_eps_min = 0.001;
  g_eps_max = 0.20;
  g_eps_values = linspace(g_eps_min,g_eps_max,ng)';
  p_g_eps_values = ones(ng,1)/ng;
  % Convert to standard deviation per period
  g_dtau_values = g_dtau_values/sqrt(nper);
  g_eps_values = g_eps_values/sqrt(nper);
  % Convert to g-values for variances instead of standard deviations (ln(s^2) = 2*ln(s))
  g_dtau_values=2*g_dtau_values;
  g_eps_values=2*g_eps_values;
  % Save priors
  g_eps_prior = [g_eps_values p_g_eps_values];
  g_dtau_prior = [g_dtau_values p_g_dtau_values];
   
  for iseries = 1:size(dp,2)
    
      str_name = char(namevec(iseries));
    ulabel = [str_name '_poos_' outlabel]; 
    fprintf(['Carrying out Calculations for ' ulabel '\n']);
    datetime('now');
    y = dp(ismpl==1,iseries);

    % Scale y:  This scale matters because of the off-set parameters "c" in 
    % draw_lcs_indicators and draw_sigma;
    dy = y(2:end)-y(1:end-1);
    inan = isnan(dy);
    std_dy = std(dy(inan==0));
    scale_y = std_dy/5;
    yn = y/scale_y;
    
    % Matrices for saving results
    tau_draws = NaN*zeros(dnobs_ismpl,n_draws_save);
    tic;
    for t = i_first_poos:dnobs_ismpl;
      fprintf([str_name '  %4i : %2i '],calds_ismpl(t,:));
      toc
      tic;
      sigma_eps = ones(t,1);
      sigma_dtau = ones(t,1);
      scale_eps = ones(t,1);
      ps = ps_prior(1)/(ps_prior(1)+ps_prior(2));
      
      % Burnin Draws
      [tau,tau_f,sigma_dtau,sigma_eps,g_eps,g_dtau,scale_eps,ps] = ucsv_outlier_draw_const_vol_e(yn(1:t),n_burnin,sigma_eps,sigma_dtau,scale_eps,ps,g_eps_prior,g_dtau_prior,scl_eps_vec,ps_prior);
      
      % MCMC Draws to save at date t
      for i = 1:n_draws_save;
       [tau,tau_f,sigma_dtau,sigma_eps,g_eps,g_dtau,scale_eps,ps] = ucsv_outlier_draw_const_vol_e(yn(1:t),k_draws,sigma_eps,sigma_dtau,scale_eps,ps,g_eps_prior,g_dtau_prior,scl_eps_vec,ps_prior);
       tau_draws(t,i) = tau(end);
      end;
    end;
     
    % Rescale so that things are in units of y
    tau_draws = tau_draws*scale_y;  
    % Compute posterior for variances
      pctvec = [0.05 1/6 0.50 5/6 0.95]';
      tau_mean_pct = post_mean_pct(tau_draws',pctvec);
      str_tmp = [matdir 'tau_mean_pct' ulabel]; save(str_tmp,'tau_mean_pct');    
  end;