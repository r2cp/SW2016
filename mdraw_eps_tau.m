function [eps_common,tau_a_common,tau_a_unique,tau_f_common,tau_f_unique] = mdraw_eps_tau(y,alpha_eps,alpha_tau,sigma_eps_common,sigma_dtau_common,sigma_eps_unique,sigma_dtau_unique);
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
  X_Draw_Filt = NaN*zeros(nobs+1,ns);  % Draw from filtered .. these are marginal draws, not joint (same as ucsv) ;
  
  % Check for missing values
  if sum(sum(isnan(y))) == 0;
  for t = 1:nobs;
    Ht = H;
    Ht(1,:) = alpha_eps(t,:);
    Ht(2,:) = alpha_tau(t,:);
    Q_t = diag([(sigma_eps_common(t,1).^2) (sigma_dtau_common(t,1).^2) (sigma_dtau_unique(t,:).^2)]');
    Rt = diag((sigma_eps_unique(t,:).^2)');
    yt = y(t,:)';
    [X1,P1,X2,P2] = kfilt(yt,X1,P1,Ht,F,Rt,Q_t);
    X1t(:,t+1)=X1;
    P1t(:,:,t+1)=P1;
    X2t(:,t+1)=X2;
    P2t(:,:,t+1)=P2;
    chol_P1 = chol(P1);
    X = X1+chol_P1'*randn(ns,1);
    X_Draw_Filt(t+1,:)=X';
   end;
  
  else;
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
    [X1,P1,X2,P2] = kfilt(yt,X1,P1,Ht,F,Rt,Q_t);
    X1t(:,t+1)=X1;
    P1t(:,:,t+1)=P1;
    X2t(:,t+1)=X2;
    P2t(:,:,t+1)=P2;
    chol_P1 = chol(P1);
    X = X1+chol_P1'*randn(ns,1);
    X_Draw_Filt(t+1,:)=X';
   end;
 
 end; 
   
   % Draw From State
   X_Draw = NaN*zeros(nobs+1,ns);
   % Initial Draw
   P3 = P1;
   X3 = X1;
   chol_P3 = chol(P3);
   X = X3+chol_P3'*randn(ns,1);
   X_Draw(nobs+1,:)=X';
   for t = nobs:-1:1;
      X1 = X1t(:,t);
      X2 = X2t(:,t+1);
      P1 = P1t(:,:,t);
      P2 = P2t(:,:,t+1);
      F_P1 = F*P1;
      [AS,rcond]=linsolve(P2,F_P1);
      if rcond < 1.0e-12;;
           P2i=pinv(P2);
           AS = P2i*F_P1;
      end;
      P3=P1-AS'*F_P1;
      P3 = 0.5*(P3+P3');
      X3=X1+AS'*(X-X2);
      X = X3;
      if t > 1;
          chol_P3 = chol(P3);    
          X = X + chol_P3'*randn(ns,1);
      else;
          P3 = P3(3:end,3:end);
          chol_P3 = chol(P3);    
          X(3:end) = X(3:end) + chol_P3'*randn(ns-2,1);
      end;             
      X_Draw(t,:)=X';
   end;
   eps_common = X_Draw(2:end,1);
   tau_a_common = X_Draw(:,2);
   tau_a_unique = X_Draw(:,3:ns);
   tau_f_common = X_Draw_Filt(2:end,2);
   tau_f_unique = X_Draw_Filt(2:end,3:ns);

end