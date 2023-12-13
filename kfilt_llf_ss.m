function [X1,X2,e,llf] = kfilt_llf_ss(y,X1,H,F,K,hti,log_det_ht);

%{
	Update of kfilt_llf.gss (GAUSS program), MWW, 10/26/2014
	Steady State Calculations
  
  Kalman Filter Procedure  -- Hamilton Notation
  
  y(t) = H'*x(t) + w(t)
  x(t) = F x(t-1) + v(t)

  var(w(t)) = R;
  var(v(t)) = Q;

  X1 = x(t-1/t-1)  -- on input
  P1 = p(t-1/t-1)  -- on input
  
  Output includes innovations, variance of innovations and log-likelihood;
%}
 
 X2=F*X1;
 e=y-H'*X2;
 X1=X2+K*e;
 llf=-0.5*(log_det_ht)+(e'*hti*e));

end