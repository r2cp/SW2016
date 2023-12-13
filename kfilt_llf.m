function [X1,P1,X2,P2,e,ht,llf] = kfilt_llf(y,X1,P1,H,F,R,Q);

%{
	Update of kfilt_llf.gss (GAUSS program), MWW, 10/26/2014
  
  Kalman Filter Procedure  -- Hamilton Notation
  
  y(t) = H'*x(t) + w(t)
  x(t) = F x(t-1) + v(t)

  var(w(t)) = R;
  var(v(t)) = Q;

  X1 = x(t-1/t-1)  -- on input
  P1 = p(t-1/t-1)  -- on input
  
  Output includes innovations, variance of innovations and log-likelihood;
%}
 nstate = size(X1,1);
 eye_ny = eye(size(y,1));
 X2=F*X1;
 e=y-H'*X2;
 P2=F*P1*F'+ Q;
 ht=H'*P2*H + R; 
 opts.POSDEF=true;
 [hti,rcond]=linsolve(ht,eye_ny);
 if rcond < 1.0e-12;
   hti = pinv(ht);
 end;
 %hti=inv(ht);
 K=P2*H*hti;
 X1=X2+K*e;
 P1=(eye(nstate)-K*H')*P2;
 P1=0.5*(P1+P1');
 llf=-0.5*(log(det(ht))+(e'*hti*e));

end