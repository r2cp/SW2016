% ucsv_agg.m -- aggregate inflation
% 9/11/2015, MWW
%


clear all;
small = 1.0e-10;
big = 1.0e+6;                  
%rng(663436);
rng(99837);
 
  % -- File Directories  
  outdir = '/Users/mwatson/Dropbox/p_comp/revision/matlab/out/';
  figdir = '/Users/mwatson/Dropbox/p_comp/revision/matlab/fig/';
  matdir = '/Users/mwatson/Dropbox/p_comp/revision/matlab/mat/';
% 
%   
  % -- Read in Data --- 
  load_data = 1;  % 1 if reloading data from Excel, etc 
  mtoq_agg = 3;   % Temporal aggregation indicator of monthly to quarterly data
  pcomp_data_calendar_m_and_q;
  % Data Series Used
  dp_agg = dp_agg_q;
  dp_agg_xfe = dp_agg_xfe_q;
  dp_agg_xe = dp_agg_xe_q;
  dp_disagg = dp_disagg_q;
  calvec = calvec_q;
  dnobs = dnobs_q;
  calds = calds_q;
  nper = 4;
  outlabel = 'Q3_run2';
   
  labvec = {'Headline Inflation';'Core XFE';'Core XE'};
  namevec = {'dp_agg';'dp_xfe';'dp_xe'};
  dp = [dp_agg dp_agg_xfe dp_agg_xe];

  % Dates
  first_date = [1960 1];
  last_date = calds(end,:);
  
  ismpl = smpl(calvec,first_date,last_date,nper);
  calvec_ismpl = calvec(ismpl==1);
  dnobs_ismpl = size(calvec_ismpl,1);
  str_tmp = [matdir 'calvec_ismpl_ucsv']; save(str_tmp,'calvec_ismpl');
  
  % Parameters for UCSV Draws
  n_burnin = 10000;     % Discarded Draws
  n_draws_save = 5000;  % Number of Draws to Save
  k_draws = 10;         % Save results every k_draws
  n_draws =  n_draws_save*k_draws;   % Total Number of Draws after burnin
  
  % Parameters for scale mixture of epsilon component
  scl_eps_vec = [1; linspace(2.0,10.0,9)'];
  ps_mean = 1 - 1/(4*nper);              % Outlier every 4 years
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
   
  for iseries = 1:1;
    tmp = char(namevec(iseries));
    ulabel = [tmp '_' outlabel]; 
    fprintf(['Carrying out Calculations for ' ulabel '\n']);
    
    y = dp(ismpl==1,iseries);
    [tau_draws,tau_f_draws,sigma_dtau_draws,sigma_eps_draws,g_eps_draws,g_dtau_draws,scl_eps_draws,sigmatotal_eps_draws,ps_draws] = ucsv_outlier(y,n_burnin,n_draws,k_draws,g_eps_prior,g_dtau_prior,scl_eps_vec,ps_prior);
    
    % Compute posteriors for trend values
    pctvec = [0.05 1/6 0.50 5/6 0.95]';
    tau_mean_pct = post_mean_pct(tau_draws',pctvec);
    str_tmp = [matdir 'tau_mean_pct' ulabel]; save(str_tmp,'tau_mean_pct');
  
  end;