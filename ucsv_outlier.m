function [tau_draws,tau_f_draws,sigma_dtau_draws,sigma_eps_draws,g_eps_draws,g_dtau_draws,scl_eps_draws,sigmatotal_eps_draws,ps_draws] = ucsv_outlier(y,n_burnin,n_draws,k_draws,g_eps_prior,g_dtau_prior,scl_eps_vec,ps_prior)

% Model and notation
  %  y(t) = tau(t) + eps(t)
  %    tau(t) = tau(t-1) + dtau(t)
  %    eps and d_tau(t) are white noises with
  %     time varying standard deviations sigma_eps(t) and sigma_dtau(t)
  %     The log-volatilities follow random walk process with innovation standard deviations given by the parameter "g"
  %     The prior for g is given in the input matrices g_xxx_prior.  The first column are a grid of values for g and the 
  %     the second column are the probabilities
  %     This version of the program allows for fat tails in the distribution of eps by using a scale mixture of normals
  %      The scales (standard deviations) are given in scl_eps_vec and their prior probabilities are given in prob_scl_eps_vec
  %
  %     Note: to impose a particular value of g, use a point-mass prior  
  %     The output below are the draws of these variables
  %      tau_draws are the draws of tau
  %      tau_f_draws are "filtered" draws of tau -- but using the full-
  %        sample draws of sigma_eps and sigma_dtau
  %      sigma_eps_draws are the draws of sigma_eps
  %      sigma_dtau_draws are the draws of sigma_dtau
  %      scl_eps_draw are the draws of the scale of epsilon (normalized by volatility)
  
  
  n_draws_save = floor(n_draws/k_draws);  % Total of n_draws, every k_draws is savee
  
  
  % Mixture parameters for log chisquared (1) approximation 
  % r_p are probabilities
  % r_m are means
  % r_s are standard deviations
  
  %{
  % 2-component mixture used in Stock-Watson (2007)
  r_p = [0.086 (1-0.086)]';
  r_m = [-7.472 -0.698]';
  r_s = [1.411 1.411]';
  %}
  
  %{
  % 7-component mixture from Kim, Shephard, and Chib ReStud (1998);
  r_p = [0.00730 0.10556 0.00002 0.04395 0.34001 0.24566 0.25750]';
  r_m = [-10.12999 -3.97281 -8.56686 2.77786 0.61942 1.79518 -1.08819]';
  r_v = [5.79596 2.61369 5.17950 0.16735 0.64009 0.34023 1.26261]';
  r_m = r_m - 1.2704;
  r_s = sqrt(r_v);
  %}
  
  % 10-component mixture from Omori, Chib, Shephard, and Nakajima JOE (2007)
  r_p = [0.00609 0.04775 0.13057 0.20674 0.22715 0.18842 0.12047 0.05591 0.01575 0.00115]';
  r_m = [1.92677 1.34744 0.73504 0.02266 -0.85173 -1.97278 -3.46788 -5.55246 -8.68384 -14.65000]';
  r_v = [0.11265 0.17788 0.26768 0.40611 0.62699 0.98583 1.57469 2.54498 4.16591 7.33342]'; 
  r_s = sqrt(r_v);
  
  T = size(y,1);   % Sample size
  
  % Scale y:  This scale matters because of the off-set parameters "c" in 
  % draw_lcs_indicators and draw_sigma;
  dy = y(2:end)-y(1:end-1);
  inan = isnan(dy);
  std_dy = std(dy(inan==0));
  scale_y = std_dy/5;
  yn = y/scale_y;
  
  % Matrices for saving results
  tau_draws = NaN*zeros(T,n_draws_save);
  tau_f_draws = NaN*zeros(T,n_draws_save);
  sigma_dtau_draws = NaN*zeros(T,n_draws_save);
  sigma_eps_draws = NaN*zeros(T,n_draws_save);
  g_eps_draws = NaN*zeros(n_draws_save,1);
  g_dtau_draws = NaN*zeros(n_draws_save,1);
  scl_eps_draws = NaN*zeros(T,n_draws_save);
  sigmatotal_eps_draws = NaN*zeros(T,n_draws_save);
  ps_draws = NaN*zeros(n_draws_save,1);
  
  % -- Begin Calculations: sigma = 1 and scale_eps is unity;
  sigma_eps = ones(T,1);
  sigma_dtau = ones(T,1);
  scale_eps = ones(T,1);
  
  % Initial value of ps is the mean
  n_scl_eps = length(scl_eps_vec);
  ps = ps_prior(1)/(ps_prior(1)+ps_prior(2));
  ps2 = (1-ps)/(n_scl_eps-1);
  ps2 = ps2*ones(n_scl_eps-1,1);
  prob_scl_eps_vec = [ps ; ps2];
  
  itmp = 0;
  tic;
  kk_draw  = 0;
  jj_draw = 0;
  
  for idraw=1:n_burnin+n_draws;
  	% Step 1(a): draw tau, tau_f, dtau, and eps; 
  	sigma_eps_scl = sigma_eps.*scale_eps;	  % SD of eps, which is stochastic volatility times scale in mixture distribution 
    [tau_a,tau_f] = draw_tau(yn,sigma_dtau,sigma_eps_scl);
    dtau = tau_a(2:end)-tau_a(1:end-1);
    tau = tau_a(2:end);
    eps = yn - tau;
    eps_scaled = eps./scale_eps;  % Scaled version of epsilon
  
    % Step 1(b): Draw mixture indicators for log chi-squared(1)
    ind_eps = draw_lcs_indicators(eps_scaled,sigma_eps,r_p,r_m,r_s); 
    ind_dtau = draw_lcs_indicators(dtau,sigma_dtau,r_p,r_m,r_s);
    
    % Step 2(a): Draw G
    g_eps = draw_g(eps_scaled,g_eps_prior,ind_eps,r_m,r_s,1);
    g_dtau = draw_g(dtau,g_dtau_prior,ind_dtau,r_m,r_s,1);
    
    % Step 2(b): Draw Volatilities
    sigma_eps = draw_sigma(eps_scaled,g_eps,ind_eps,r_m,r_s,1);	
    sigma_dtau = draw_sigma(dtau,g_dtau,ind_dtau,r_m,r_s,1);	
    
    % Step 3: Draw Scale of epsilon
    % scale_eps = draw_scale_eps(eps,sigma_eps,ind_eps,r_m,r_s,scl_eps_vec,prob_scl_eps_vec);
    % Scale for epsilon = 1
    scale_eps = ones(T, 1); 

    % Step 4; Draw probability of outlier;
    % prob_scl_eps_vec = draw_ps(scale_eps,ps_prior,n_scl_eps);
    % Remove probability of outlier
    prob_scl_eps_vec(1) = 1; 
    prob_scl_eps_vec(2:end) = 0; 
   
    % Save draws;
    if idraw > n_burnin;
      kk_draw = kk_draw+1;
      if kk_draw == k_draws;
          jj_draw = jj_draw+1;
          kk_draw = 0;
          tau_draws(:,jj_draw) = tau;
          tau_f_draws(:,jj_draw) = tau_f;
          sigma_dtau_draws(:,jj_draw) = sigma_dtau;
          sigma_eps_draws(:,jj_draw) = sigma_eps;
          g_eps_draws(jj_draw) = g_eps;
          g_dtau_draws(jj_draw) = g_dtau;
          scl_eps_draws(:,jj_draw) = scale_eps;
          sigmatotal_eps_draws(:,jj_draw) = sigma_eps.*scale_eps;
          ps_draws(jj_draw) = prob_scl_eps_vec(1);
      end;
    end;
    
    itmp=itmp+1;
    if itmp == 10000;
      itmp=0;
      fprintf('  Number of Draws: %6i    ',idraw);
      toc
      tic;
     end;
    
  end;
  
  % Rescale so that things are in units of y
  tau_draws = tau_draws*scale_y;
  tau_f_draws = tau_f_draws*scale_y;
  sigma_dtau_draws = sigma_dtau_draws*scale_y;
  sigma_eps_draws = sigma_eps_draws*scale_y;
  sigmatotal_eps_draws = sigmatotal_eps_draws*scale_y;

end

