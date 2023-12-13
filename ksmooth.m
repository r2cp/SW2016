function [X3,P3] = ksmooth(X1,X2,X3,P1,P2,P3,F);
 tmp = F*P1;
 [AS,R]=linsolve(P2,tmp);
 opts.POSDEF=true;
 if R < 1.0e-12;;
  P2i=pinv(P2);
  AS = P2i*tmp;
 end;
 P3=P1+AS'*(P3-P2)*AS;
 P3 = 0.5*(P3+P3');
 X3=X1+AS'*(X3-X2);
end
