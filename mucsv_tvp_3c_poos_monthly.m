% mucsv_tvp.m
% Updated Sept 28, 2015
% pcomp data
% multivariate UCSV Model 
% TVP (Random Walk) Factor loadings

clear all;
small = 1.0e-10;
big = 1.0e+6;                  
rng(663436);

  % -- File Directories  
  outdir = '/Users/mwatson/Dropbox/p_comp/revision/matlab/out/';
  figdir = '/Users/mwatson/Dropbox/p_comp/revision/matlab/fig/';
  matdir = '/Users/mwatson/Dropbox/p_comp/revision/matlab/mat/';
  
 % -- Read in Data --- 
  load_data = 1;  % 1 if reloading data from Excel, etc 
  mtoq_agg = 0;   % Temporal aggregation indicator of monthly to quarterly data
  pcomp_data_calendar_m_and_q;  
  % Data Series Used
  dp_agg = dp_agg_m;
  dp_agg_xfe = dp_agg_xfe_m;
  dp_agg_xe = dp_agg_xe_m;
  dp_disagg = dp_disagg_3comp_m;
  share_avg = share_avg_3comp_m;
  share_avg_xfe = share_avg_3comp_xfe_m;
  share_avg_xe = share_avg_3comp_xe_m;
  calvec = calvec_m;
  dnobs = dnobs_m;
  calds = calds_m;
  nper = 12;
  mlabel = '_mucsv_tvp_3c_M';
  
  labelvec_disagg = labelvec_disagg_3comp;
  namevec_disagg = namevec_disagg_3comp;
  n_incl = n_incl_3comp;
  
  % --- Sample Period for Analyais 
  first_date = [1960 1];
  last_date = calds(end,:);
  ismpl = smpl(calvec,first_date,last_date,nper);
  calvec_ismpl = calvec(ismpl==1);
  calds_ismpl = calds(ismpl==1,:);
  notim = size(calvec_ismpl,1);
  dnobs_ismpl = notim;
  str_tmp = [matdir 'calvec_ismpl' mlabel]; save(str_tmp,'calvec_ismpl');
  
  % Chain Weights over this sample period 
  cw = share_avg(ismpl==1,:);
  cw_xfe = share_avg_xfe(ismpl==1,:);
  cw_xe = share_avg_xe(ismpl==1,:);
  
   % First POOS Date
  first_date_poos = [1990 1];
  tmp = calds_ismpl == repmat(first_date_poos,dnobs_ismpl,1);
  [a,i_first_poos] = max(sum(tmp,2));  % Index in calvec_ismpl of first poos period
 
  % Parameters for UCSV Draws -- First Pass
  n_burnin = 2000;      % Discarded Draws
  n_draws_save = 2000;  % Number of Draws to Save
  k_draws = 5;          % Save results every k_draws
 
  % Parameters for scale mixture of epsilon component
  scl_eps_vec = [1; linspace(2.0,10.0,9)'];
  ps_mean = 1 - 1/(4*nper);              % Outlier every 4 years
  ps_prior_obs = nper*10;                % Sample size of 10 years for prior
  ps_prior_a = ps_mean*ps_prior_obs;     % "alpha" in beta prior
  ps_prior_b = (1-ps_mean)*ps_prior_obs; % "beta" in beta prior
  ps_prior = [ps_prior_a ps_prior_b];
  
  % -- Data -- 
  dp_mat = dp_disagg(ismpl==1,:);
  n_y = size(dp_mat,2);
   
  % Scale data:  This scale matters because of the off-set parameters "c" in 
  % draw_lcs_indicators and draw_sigma;
  dy = dp_mat(2:end,:)-dp_mat(1:end-1,:);
  sd_ddp = NaN*zeros(n_y,1);
  for i = 1:n_y;
   tmp = dy(:,i);
   inan = isnan(tmp);
   sd_ddp(i) = std(tmp(inan==0));
  end;
  sd_ddp_median = median(sd_ddp);
  scale_y = sd_ddp_median/5;
  dp_mat_n = dp_mat/scale_y;
  
  % -- Parameters for model
   % 10-component mixture approximation to log chi-squared(1) from Omori, Chib, Shephard, and Nakajima JOE (2007)
   r_p = [0.00609 0.04775 0.13057 0.20674 0.22715 0.18842 0.12047 0.05591 0.01575 0.00115]';
   r_m = [1.92677 1.34744 0.73504 0.02266 -0.85173 -1.97278 -3.46788 -5.55246 -8.68384 -14.65000]';
   r_v = [0.11265 0.17788 0.26768 0.40611 0.62699 0.98583 1.57469 2.54498 4.16591 7.33342]'; 
   r_s = sqrt(r_v);
   
   % Prior for g
   % Gvalues .. for evolution of standard deviations over a year (nper periods)
   ng = 5;      % Number of grid points for approximate uniform prior
   g_dtau_min = 0.001;
   g_dtau_max = 0.20;
   g_dtau_unique_values = linspace(g_dtau_min,g_dtau_max,ng)';
   p_g_dtau_unique_values = ones(ng,1)/ng;
   
   g_eps_min = 0.001;
   g_eps_max = 0.20;
   g_eps_unique_values = linspace(g_eps_min,g_eps_max,ng)';
   p_g_eps_unique_values = ones(ng,1)/ng;
   
   g_dtau_min = 0.001;
   g_dtau_max = 0.20;
   g_dtau_common_values = linspace(g_dtau_min,g_dtau_max,ng)';
   p_g_dtau_common_values = ones(ng,1)/ng;
   
   g_eps_min = 0.001;
   g_eps_max = 0.20;
   g_eps_common_values = linspace(g_eps_min,g_eps_max,ng)';
   p_g_eps_common_values = ones(ng,1)/ng;
   
   
   % Convert to standard deviation per period
   g_dtau_unique_values = g_dtau_unique_values/sqrt(nper);
   g_eps_unique_values = g_eps_unique_values/sqrt(nper);
   g_dtau_common_values = g_dtau_common_values/sqrt(nper);
   g_eps_common_values = g_eps_common_values/sqrt(nper);
   
   % Convert to g-values for variances instead of standard deviations (ln(s^2) = 2*ln(s))
   g_dtau_unique_values=2*g_dtau_unique_values;
   g_eps_unique_values=2*g_eps_unique_values;
   g_dtau_common_values=2*g_dtau_common_values;
   g_eps_common_values=2*g_eps_common_values;
   
   % Save priors
   g_eps_unique_prior = [g_eps_unique_values p_g_eps_unique_values];
   g_dtau_unique_prior = [g_dtau_unique_values p_g_dtau_unique_values];
   g_eps_common_prior = [g_eps_common_values p_g_eps_common_values];
   g_dtau_common_prior = [g_dtau_common_values p_g_dtau_common_values];
   
   % Parameters for prior for factor loadings -- note these
   % depend on scaling (scale_y) introduced above, so they are in scaled y
   % units
   % ... inital values of factor loadings
   omega_tau = 10/scale_y;
   omega_eps = 10/scale_y;
   sigma_tau = 0.4/scale_y;
   sigma_eps = 0.4/scale_y;
   var_alpha_tau = ((omega_tau^2)*ones(n_y,n_y)) + ((sigma_tau^2)*eye(n_y));
   var_alpha_eps = ((omega_eps^2)*ones(n_y,n_y)) + ((sigma_eps^2)*eye(n_y));
   prior_var_alpha = zeros(2*n_y,2*n_y);
   prior_var_alpha(1:n_y,1:n_y) = var_alpha_eps;
   prior_var_alpha(n_y+1:2*n_y,n_y+1:2*n_y) = var_alpha_tau;
   % Alpha TVP parameters -- use "Number of prior obs (nu) and prior squared (s2)" as parameters .. as in Del Negro and Otrok;
   nu_prior_alpha = 0.1*notim;             
   s2_prior_alpha = (0.25/sqrt(notim))^2;
   s2_prior_alpha = s2_prior_alpha/(scale_y^2);
   
   
   % Matrix for saving results to compute final values
   agg_tau_draws = NaN*ones(notim,n_draws_save);
   agg_tau_xfe_draws = NaN*ones(notim,n_draws_save);
   agg_tau_xe_draws = NaN*ones(notim,n_draws_save);
   
   tic;
   for t = i_first_poos:notim;
    fprintf('%4i : %2i ',calds_ismpl(t,:));
    toc
    tic;
   
    % Initial Values of parameters -- I save alpha values for each date, so the
    % same program can be used when TVP is allowed
    %alpha_eps = rand(t,n_y);
    %alpha_tau = rand(t,n_y);
    alpha_eps = ones(t,n_y);
    alpha_tau = ones(t,n_y);
    sigma_dtau_unique = rand(t,n_y);
    sigma_eps_unique = rand(t,n_y);
    sigma_dtau_common = rand(t,1);
    sigma_eps_common = rand(t,1);
    scale_eps_unique = ones(t,n_y);
    scale_eps_common = ones(t,1);
    % Initial values for sigma_alpha
    a = nu_prior_alpha/2;
    ssr = nu_prior_alpha*s2_prior_alpha;
    b = 2/ssr;
    var_dalpha = 1./gamrnd(a,b,[2*n_y,1]);
    sigma_dalpha = sqrt(var_dalpha);
    % Initial Value of ps
    n_scl_eps = length(scl_eps_vec);
    ps = ps_prior(1)/(ps_prior(1)+ps_prior(2));
    ps2 = (1-ps)/(n_scl_eps-1);
    ps2 = ps2*ones(n_scl_eps-1,1);
    prob_scl_eps_vec = [ps ; ps2];
    prob_scl_eps_vec_common = prob_scl_eps_vec;
    prob_scl_eps_vec_unique = repmat(prob_scl_eps_vec,1,n_y);
   
    for jj_draw = 0:n_draws_save;
     if jj_draw == 0;
       kk = n_burnin;
     else;
       kk = k_draws;
     end;
     for i_draw = 1:kk;
  	    sigma_eps_unique_scl = sigma_eps_unique.*scale_eps_unique;  % SD of eps_unique, which is stochastic volatility times scale in mixture distribution 
  	    sigma_eps_common_scl = sigma_eps_common.*scale_eps_common;  % SD of eps_common, which is stochastic volatility times scale in mixture distribution 
  	    % Step 1.a.1: draw tau, tau_f, dtau, and eps; 	
        [eps_common,tau_a_common,tau_a_unique,tau_f_common,tau_f_unique] = mdraw_eps_tau(dp_mat_n(1:t,:),alpha_eps,alpha_tau,sigma_eps_common_scl,sigma_dtau_common,sigma_eps_unique_scl,sigma_dtau_unique);
        dtau_common = tau_a_common(2:end)-tau_a_common(1:end-1);
        dtau_unique = tau_a_unique(2:end,:)-tau_a_unique(1:end-1,:);
        tau_common = tau_a_common(2:end);
        tau_unique = tau_a_unique(2:end,:);
    
        % Step 1.a.2 : Draw Factor Loadings 
        % -- Draw alpha_eps and alpha_tau;
        [alpha_eps, alpha_tau, dalpha_eps, dalpha_tau] = draw_alpha_tvp(dp_mat_n(1:t,:),prior_var_alpha,sigma_dalpha,tau_unique,eps_common,tau_common,sigma_eps_unique_scl); 
     
        % Step 1.a.3: Draw Standard Deviations of Alpha TVPs;
        dalpha = [dalpha_eps dalpha_tau];
        sigma_dalpha = draw_dalpha_sigma(dalpha,nu_prior_alpha, s2_prior_alpha);
    
        % Save some values
        y_eps_common = alpha_eps.*repmat(eps_common,1,n_y);
        y_tau_common = alpha_tau.*repmat(tau_common,1,n_y);
        y_tau_unique = tau_unique;
        eps_unique = dp_mat_n(1:t,:) - y_eps_common - y_tau_common-y_tau_unique;
        y_eps_unique = eps_unique;
     
        % Step 1(b): Draw mixture indicators for log chi-squared(1)
        eps_unique = dp_mat_n(1:t,:) - y_eps_common - y_tau_common-y_tau_unique;
        eps_unique_scaled = eps_unique./scale_eps_unique;
        eps_common_scaled = eps_common./scale_eps_common;
        ind_eps_common = draw_lcs_indicators(eps_common_scaled,sigma_eps_common,r_p,r_m,r_s); 
        ind_dtau_common = draw_lcs_indicators(dtau_common,sigma_dtau_common,r_p,r_m,r_s);
        ind_eps_unique = NaN*ones(t,size(r_p,1),n_y);
        ind_dtau_unique = NaN*ones(t,size(r_p,1),n_y);
        for i = 1:n_y;
         tmp = draw_lcs_indicators(eps_unique_scaled(:,i),sigma_eps_unique(:,i),r_p,r_m,r_s); 
         ind_eps_unique(:,:,i) = tmp;
         tmp = draw_lcs_indicators(dtau_unique(:,i),sigma_dtau_unique(:,i),r_p,r_m,r_s); 
         ind_dtau_unique(:,:,i) = tmp;
        end;
     
        % Step 2(a): Draw G
        i_init = 0;   % Variance = 1 as initial condition to identify factor loadings
        g_eps_common = draw_g(eps_common_scaled,g_eps_common_prior,ind_eps_common,r_m,r_s,i_init);
        g_dtau_common = draw_g(dtau_common,g_dtau_common_prior,ind_dtau_common,r_m,r_s,i_init);
    
        i_init = 1;   % Vague prior for initial variance;
        g_eps_unique = NaN*zeros(n_y,1);
        g_dtau_unique = NaN*zeros(n_y,1);
        for i = 1:n_y;
          g_eps_unique(i) = draw_g(eps_unique_scaled(:,i),g_eps_unique_prior,ind_eps_unique(:,:,i),r_m,r_s,i_init);
          g_dtau_unique(i) = draw_g(dtau_unique(:,i),g_dtau_unique_prior,ind_dtau_unique(:,:,i),r_m,r_s,i_init);
        end;
     
        % Step 2(b): Draw Volatilities
        i_init = 0;   % Variance = 1 as initial condition to identify factor loadings
        sigma_eps_common = draw_sigma(eps_common_scaled,g_eps_common,ind_eps_common,r_m,r_s,i_init);	
        sigma_dtau_common = draw_sigma(dtau_common,g_dtau_common,ind_dtau_common,r_m,r_s,i_init);	
    
        i_init = 1;  % Vague prior for initial variance;
        sigma_eps_unique = NaN*ones(t,n_y);
        sigma_dtau_unique = NaN*ones(t,n_y);
        for i = 1:n_y;
          sigma_eps_unique(:,i) = draw_sigma(eps_unique_scaled(:,i),g_eps_unique(i),ind_eps_unique(:,:,i),r_m,r_s,i_init);  
          sigma_dtau_unique(:,i) = draw_sigma(dtau_unique(:,i),g_dtau_unique(i),ind_dtau_unique(:,:,i),r_m,r_s,i_init);
        end;
    
        % Step 3: Draw Scale of epsilon
        scale_eps_common = draw_scale_eps(eps_common,sigma_eps_common,ind_eps_common,r_m,r_s,scl_eps_vec,prob_scl_eps_vec_common);
        scale_eps_unique = NaN*ones(t,n_y);
        for i = 1:n_y; 
          scale_eps_unique(:,i) = draw_scale_eps(eps_unique(:,i),sigma_eps_unique(:,i),ind_eps_unique(:,:,i),r_m,r_s,scl_eps_vec,prob_scl_eps_vec_unique(:,i));
        end;
   
        % Step 4; Draw probability of outlier;
        prob_scl_eps_vec_common = draw_ps(scale_eps_common,ps_prior,n_scl_eps);
        for i = 1:n_y;
         prob_scl_eps_vec_unique(:,i) = draw_ps(scale_eps_unique(:,i),ps_prior,n_scl_eps);
        end;
      end;
      if jj_draw > 0;
        agg_tau_draws(t,jj_draw) = cw(t,:)*(y_tau_common(t,:)+y_tau_unique(t,:))';
        agg_tau_xfe_draws(t,jj_draw) = cw_xfe(t,:)*(y_tau_common(t,:)+y_tau_unique(t,:))';
        agg_tau_xe_draws(t,jj_draw) = cw_xe(t,:)*(y_tau_common(t,:)+y_tau_unique(t,:))';
      end;
   end;
 end;
     
 agg_tau_draws = agg_tau_draws*scale_y;   % Rescale
 agg_tau_xfe_draws = agg_tau_xfe_draws*scale_y;   % Rescale
 agg_tau_xe_draws = agg_tau_xe_draws*scale_y;   % Rescale
  
 pctvec = [0.05 1/6 0.50 5/6 0.95]';
 
 agg_tau_mean_pct = post_mean_pct(agg_tau_draws',pctvec);
 str_tmp = [matdir 'agg_tau_mean_pct_poos' mlabel]; save(str_tmp,'agg_tau_mean_pct');
 
 agg_tau_xfe_mean_pct = post_mean_pct(agg_tau_xfe_draws',pctvec);
 str_tmp = [matdir 'agg_tau_xfe_mean_pct_poos' mlabel]; save(str_tmp,'agg_tau_xfe_mean_pct');
 
 agg_tau_xe_mean_pct = post_mean_pct(agg_tau_xe_draws',pctvec);
 str_tmp = [matdir 'agg_tau_xe_mean_pct_poos' mlabel]; save(str_tmp,'agg_tau_xe_mean_pct');
  