% ucsv_agg.m -- aggregate inflation
% 9/11/2015, MWW
%


clear all;
small = 1.0e-10;
big = 1.0e+6;                  
rng(663436);
 
  % -- File Directories  
  outdir = 'matlab/out/';
  figdir = 'matlab/fig/';
  matdir = 'matlab/mat/mucsv_tvp_17c_Q3/';
%     
  
  % -- Read in Data --- 
  load_data = 1;  % 1 if reloading data from Excel, etc 
  mtoq_agg = 3;   % Temporal aggregation indicator of monthly to quarterly data
  
  % Load data script
  pcomp_data_calendar_m_and_q;

  %% Stating data series used in this script 
  % Data Series Used

  % Aggregate series 
  dp_agg = dp_agg_q; % Quarterly

  % Excluding food and energy (fe) 
  dp_agg_xfe = dp_agg_xfe_q; % Quarterly 
  
  % Excluding energy (e) 
  dp_agg_xe = dp_agg_xe_q; % Quarterly 

  % This I don't understand 
  % dp_disagg = dp_disagg_q; % Unused in this script 
  calvec = calvec_q;
  dnobs = dnobs_q;
  calds = calds_q;
  nper = 4;
  outlabel = 'Q3';
   
  % Estimate the model for each of these series: 
  labvec = {'Headline Inflation';'Core XFE';'Core XE'};
  namevec = {'dp_agg';'dp_xfe';'dp_xe'};
  % dp = [dp_agg dp_agg_xfe dp_agg_xe];
  dp = [dp_agg]; % Only aggregate? 

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
   
  for iseries = 1:size(dp,2);
    tmp = char(namevec(iseries));
    ulabel = [tmp '_' outlabel]; 
    fprintf(['Carrying out Calculations for ' ulabel '\n']);
    
    y = dp(ismpl==1,iseries);
    [tau_draws,tau_f_draws,sigma_dtau_draws,sigma_eps_draws,g_eps_draws,g_dtau_draws,scl_eps_draws,sigmatotal_eps_draws,ps_draws] = ucsv_outlier(y,n_burnin,n_draws,k_draws,g_eps_prior,g_dtau_prior,scl_eps_vec,ps_prior);
    
    % IMA(1,1) Parameters
    lam0=(sigma_dtau_draws.^2)+2*(sigma_eps_draws.^2);
    rho1=-(sigma_eps_draws.^2)./lam0;
    theta_draws=-(ones(size(rho1))-sqrt(ones(size(rho1))-4*rho1.^2))./(2*rho1);
    var_a=lam0./(ones(size(rho1))+theta_draws.^2);
    sigma_a_draws=sqrt(var_a);
 
    % Compute Posterior for g's
    g_eps_prior_post = g_summary(g_eps_prior,g_eps_draws);
    g_dtau_prior_post = g_summary(g_dtau_prior,g_dtau_draws);
    str_tmp = [matdir 'g_eps_prior_post' ulabel]; save(str_tmp,'g_eps_prior_post');
    str_tmp = [matdir 'g_dtau_prior_post' ulabel]; save(str_tmp,'g_dtau_prior_post');
 
    % Compute posterior for variances
    pctvec = [0.05 1/6 0.50 5/6 0.95]';
    tmp = sigma_dtau_draws.^2;
    var_dtau_mean_pct = post_mean_pct(tmp',pctvec);
    tmp = sigma_eps_draws.^2;
    var_eps_mean_pct = post_mean_pct(tmp',pctvec);
    tmp = sigmatotal_eps_draws.^2;
    vartotal_eps_mean_pct = post_mean_pct(tmp',pctvec);
    str_tmp = [matdir 'var_eps_mean_pct' ulabel]; save(str_tmp,'var_eps_mean_pct');
    str_tmp = [matdir 'var_dtau_mean_pct' ulabel]; save(str_tmp,'var_dtau_mean_pct');
    str_tmp = [matdir 'vartotal_eps_mean_pct' ulabel]; save(str_tmp,'vartotal_eps_mean_pct');
 
    % Compute posteriors for trend values
    tau_mean_pct = post_mean_pct(tau_draws',pctvec);
    tau_f_mean_pct = post_mean_pct(tau_f_draws',pctvec);
    str_tmp = [matdir 'tau_mean_pct' ulabel]; save(str_tmp,'tau_mean_pct');
    str_tmp = [matdir 'tau_f_mean_pct' ulabel]; save(str_tmp,'tau_f_mean_pct');
    tau_f_mean = mean(tau_f_draws,2);
    tmp = mean(tau_f_draws.^2,2);
    tau_f_var = tmp-tau_f_mean.^2;
    tau_f_sd = sqrt(tau_f_var);
    tau_f_mean_sd = [tau_f_mean tau_f_sd];
    str_tmp = [matdir 'tau_f_mean_sd' ulabel]; save(str_tmp,'tau_f_mean_sd');
  
    tau_mean = mean(tau_draws,2);
    tmp = mean(tau_draws.^2,2);
    tau_var = tmp-tau_mean.^2;
    tau_sd = sqrt(tau_var);
    tau_mean_sd = [tau_mean tau_sd];
    str_tmp = [matdir 'tau_mean_sd' ulabel]; save(str_tmp,'tau_mean_sd');
 
    % Compute posteriors for MA parameter
    theta_mean_pct = post_mean_pct(theta_draws',pctvec);
    sigma_a_mean_pct = post_mean_pct(sigma_a_draws',pctvec);
    str_tmp = [matdir 'theta_mean_pct' ulabel]; save(str_tmp,'theta_mean_pct');
    str_tmp = [matdir 'sigma_a_mean_pct' ulabel]; save(str_tmp,'sigma_a_mean_pct');
 
    % Compute posteriors for scale_eps_values
    pctvec = [1/6 0.50 5/6]';
    scl_eps_pct = post_mean_pct(scl_eps_draws',pctvec);
    str_tmp = [matdir 'scl_eps_pct' ulabel]; save(str_tmp,'scl_eps_pct');
    
    % Compute Posterior for ps's
    ps_mean_pct = post_mean_pct(ps_draws,pctvec);
    str_tmp = [matdir 'ps_mean_pct' ulabel]; save(str_tmp,'ps_mean_pct');
    
  end;