function [tau_a,tau_f] = draw_tau(y,sigma_dtau,sigma_eps)
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
var_dtau = sigma_dtau.^2;
var_eps = sigma_eps.^2;
% P matrices 
 nobs=size(y,1);
 p1t=zeros(nobs+1,1);
 p2t=zeros(nobs+1,1);
 x1t=zeros(nobs+1,1);
 x2t=zeros(nobs+1,1);
 tau_a=zeros(nobs+1,1);
 tau_f = zeros(nobs,1);
 
% KF Using Special Structure of Problem 
 x1=0;
 p1=big;        % vague prior for initial condition 
 x1t(1,1)=x1; 
 p1t(1,1)=p1;
 
 if sum(isnan(y)) == 0;
   
    % No missing values;
    for t=1:nobs;
 	 x2=x1;
     p2=p1+var_dtau(t);
     ht=p2+var_eps(t);
     k=p2/ht;
     x1=x2+k*(y(t)-x2);
     p1=p2-k*p2;
     x1t(t+1,1)=x1;
     p1t(t+1,1)=p1;
     x2t(t+1,1)=x2;
     p2t(t+1,1)=p2;
     % Generate random draw from filtered mean and variance 
     tau_f(t) = x1 + sqrt(p1)*randn(1,1);
    end;
 
 else;
   % Check for and skip missing values 
   for t=1:nobs;
 	 x2=x1;
     p2=p1+var_dtau(t);
     if isnan(y(t)) == 1;
 	   x1 = x2;
 	   p1 = p2;
     else;
       ht=p2+var_eps(t);
       k=p2/ht;
       x1=x2+k*(y(t)-x2);
       p1=p2-k*p2;
     end;
     x1t(t+1,1)=x1;
     p1t(t+1,1)=p1;
     x2t(t+1,1)=x2;
     p2t(t+1,1)=p2;
     % Generate random draw from filtered mean and variance 
     tau_f(t) = x1 + sqrt(p1)*randn(1,1);
    end;
   
 end;
 
 % Generate Random Draws from Smoothed Distribution 
 x3mean=x1;
 p3=p1;
 chol_p=sqrt(p3);
 x3=x3mean+chol_p*randn(1,1);
 tau_a(nobs+1,1)=x3;
 for t=nobs:-1:2;
  x2=x2t(t+1);
  p2=p2t(t+1);
  x1=x1t(t);
  p1=p1t(t);
  p2i=1/p2;
  k=p1*p2i;
  e=x3-x2;
  x3mean=x1+k*e;
  p3=(1-k)*p1;
  chol_p=sqrt(p3);
  x3=x3mean+chol_p*randn(1,1);
  tau_a(t,1)=x3; 
 end;

end