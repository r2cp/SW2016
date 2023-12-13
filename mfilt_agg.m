function [agg_tau_f_total agg_tau_f_common agg_tau_f_unique] = mfilt_agg(y,alpha_eps,alpha_tau,sigma_eps_common,sigma_dtau_common,sigma_eps_unique,sigma_dtau_unique,cw);
% 
%       y: T x 1
%       sigma_dtau: T x 1
%       sigma_eps: T x 1
%       
%       Output
%       tau_a: (T+1)x1 -- tau draw .. (t = 0 through T)
%       tau_f: Tx1 -- tau_filtered draw (t = 1 through T)
      
small=1.0e-06;
big = 1.0e+6;
n_y = size(y,2);
nobs = size(y,1);

% Set up State Vector
  % --- State Vector
  %     (1) eps(t)
  %     (2) tau(t)
  %     (3) tau_u(t)
  ns = 2+n_y;      % size of state
  Q = zeros(ns,ns);
  F = zeros(ns,ns);
  F(2:ns,2:ns) = eye(ns-1);
  H = zeros(ns,n_y);
  H(3:ns,1:n_y) = eye(n_y);
  
  % Set up KF to run
  % Initial conditions
  X1_init = zeros(ns,1);
  P1_init = zeros(ns,ns);
  P1_init(3:ns,3:ns) = big*eye(n_y);  % Vague prior for tau_unique initial values 

  X1 = X1_init;
  P1 = P1_init;
  X1t=zeros(ns,nobs+1);
  P1t=zeros(ns,ns,nobs+1);
  X2t=zeros(ns,nobs+1);
  P2t=zeros(ns,ns,nobs+1);
  X1t(:,1)=X1; 
  P1t(:,:,1)=P1;
  for t = 1:nobs;
    H_t = H;
    H_t(1,:) = alpha_eps(t,:);
    H_t(2,:) = alpha_tau(t,:);
    Q_t = diag([(sigma_eps_common(t,1).^2) (sigma_dtau_common(t,1).^2) (sigma_dtau_unique(t,:).^2)]');
    R_t = diag((sigma_eps_unique(t,:).^2)');
    tmp = y(t,:)';
    inan = isnan(tmp);
    yt = tmp(inan==0);
    Ht = H_t(:,inan==0);
    Rt = R_t(inan==0,inan==0);
    [X1,P1,X2,P2,e,ht,llf] = kfilt_llf(yt,X1,P1,Ht,F,Rt,Q_t);
    X1t(:,t+1)=X1;
    P1t(:,:,t+1)=P1;
    X2t(:,t+1)=X2;
    P2t(:,:,t+1)=P2;
  end;
  
  tau_f_common = X1t(2,2:end)';
  tau_f_unique = X1t(3:ns,2:end)';
  tau_f_common = alpha_tau.*repmat(tau_f_common,1,n_y);
  tmp = tau_f_common.*cw;
  agg_tau_f_common = sum(tmp,2);
  tmp = tau_f_unique.*cw;
  agg_tau_f_unique = sum(tmp,2);
  agg_tau_f_total = agg_tau_f_common + agg_tau_f_unique;
  
  

end