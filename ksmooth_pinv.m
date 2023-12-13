function [X3,P3] = ksmooth(X1,X2,X3,P1,P2,P3,F);

 % Uses pseudo inverse for singluar P2
 P2i=pinv(P2);
 as=P1*F'*P2i;
 P3=P1+as*(P3-P2)*as';
 X3=X1+as*(X3-X2);

end
