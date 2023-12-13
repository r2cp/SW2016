function [alpha_eps, alpha_tau, dalpha_eps, dalpha_tau] = draw_alpha_tvp(y,prior_var_alpha,sigma_dalpha,tau_unique,eps_common,tau_common,sigma_eps_unique)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
     n_y = size(y,2);
     nobs = size(y,1);
     y = y - tau_unique;  % Eliminates tau_unique from y ;
     
     % Brute Force Calculation
     n_state = 2*n_y;
     % First n_y elements of state are alpha_eps; second n_y elements of state are alpha_tau
     Q = diag(sigma_dalpha.^2);
     F = eye(n_state);
     % Set up KF to run
     % Initial conditions
     X1_init = zeros(n_state,1);
     P1_init = prior_var_alpha;
     X1 = X1_init;
     P1 = P1_init;
     X1t=zeros(n_state,nobs+1);
     P1t=zeros(n_state,n_state,nobs+1);
     X2t=zeros(n_state,nobs+1);
     P2t=zeros(n_state,n_state,nobs+1);
     X1t(:,1)=X1; 
     P1t(:,:,1)=P1;
     H = zeros(n_state,n_y);
     for t = 1:nobs;
      yt = y(t,:)';
      H(1:n_y,:) = eps_common(t,1)*eye(n_y);
      H(n_y+1:2*n_y,:) = tau_common(t,1)*eye(n_y);
      Rt = diag((sigma_eps_unique(t,:).^2)');
      [X1,P1,X2,P2] = kfilt(yt,X1,P1,H,F,Rt,Q);
      X1t(:,t+1)=X1;
      P1t(:,:,t+1)=P1;
      X2t(:,t+1)=X2;
      P2t(:,:,t+1)=P2;
     end;
     % Draw From State
     X_Draw = NaN*zeros(nobs+1,n_state);
     % Initial Draw
     P3 = P1;
     X3 = X1;
     chol_P3 = chol(P3);
     X = X3+chol_P3'*randn(n_state,1);
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
      chol_P3 = chol(P3);
      X = X3+chol_P3'*randn(n_state,1);       
      X_Draw(t,:)=X';
  end;
  dalpha_eps = X_Draw(2:end,1:n_y)-X_Draw(1:end-1,1:n_y);
  dalpha_tau = X_Draw(2:end,n_y+1:2*n_y)-X_Draw(1:end-1,n_y+1:2*n_y);
  alpha_eps = X_Draw(2:end,1:n_y);
  alpha_tau = X_Draw(2:end,n_y+1:2*n_y);
end