function [betahat, seps, resid, Q, M, G] = varest(y,nlag,i_const,calvec,nfirst,nlast,nper,icomp);

% Computes VAR and covariance matrix of estimated parameters
%    
%    Input:
%    		    s = txn matrix of series
%    		    nlag = number of lags
%    		    icomp = 0 Do not compute Companion Matrices
%    		            1 Compute companion matrix ignoring contsant term
%    		            2 Compute Companion matrix adding constant as last element of state -- loadings are zero if iconst = 0
%    		    i_const = 0 Do not include constant
%    		              1 Include constant
%    		   
%                Note: icomp = 2 ... constant term is added to companion state vector (if i_const == 1)
%    		            
%    Output:
%           betahat = estimated VAR coefficients (including constant if i_const == 1)
%                  each column has coefficients for a single equation
%           seps = innovation covariance matrix
%           eps =  VAR residuals
%           nobs = number of observations used for estimation
%           Companion Form Parameter, computed in icomp = 1;
%           Q, M, G for model written in SS form
%                y(t) = Q*z(t)
%                z(t) = M*z(t-1) + G*u(t)
%                var(u(t)) = I
% 

%  Set Up VAR %

ns = size(y,2);
T = size(y,1);
x=ones(T,1);
for i = 1:nlag;
	x = [x,lag(y,i)];
end;
if i_const == 0;
 x=x(:,2:end);  % Eliminate Constant if i_const == 0;
end;

ismpl = smpl(calvec,nfirst,nlast,nper);
trnd = (1:1:T)';
trnd = trnd(ismpl==1);
y = y(ismpl==1,:);
x = x(ismpl==1,:);
tmp = packr([y x trnd]);
y = tmp(:,1:ns);
x = tmp(:,ns+1:end-1);
trnd = tmp(:,end);

betahat = x\y;
e = y - x*betahat;
ndf=size(x,1)-size(x,2);
seps=(e'*e)/ndf;
resid = NaN*zeros(T,ns);
resid(trnd,:)=e;

if icomp > 0;
   %  
   %    Transform the VAR so that it is written in Standard form as:
   %    s(t)=P1*s(t-1) + P2*s(t-2) + ... + Pvarlag*s(t-varlag) + eps(t)
   %  

   % ---- Calculate Companion Matrix ---- ;
   if i_const == 0;
      b = betahat';
      const_coef = zeros(ns,1);
   elseif i_const == 1;
      b=betahat(2:end,:)';
      const_coef = betahat(1,:)';      % Coefficients on Constant Term ;
   else;
      error('Invalid value of i_const in VAREST');
   end;
   comp=zeros(size(b,2),size(b,2));
   comp(1:size(b,1),:)=b;
   if size(b,2) > size(b,1);
     comp(size(b,1)+1:end,1:end-size(b,1))=eye(size(b,2)-size(b,1));
   end;
   %   -- Write Model in SS Form --
   %      y(t) = Q*z(t)
   %      z(t) = M*z(t-1) + G*u(t)
   %      var(u(t)) = I
   %  
   M=comp;
   Q=zeros(ns,size(M,1));
   Q(1:ns,1:ns)=eye(ns);
   G=zeros(size(M,1),ns);
   G(1:ns,1:ns)=(chol(seps))';
end;

if icomp == 2;  % Add Constant Term 
 G=[G ; zeros(1,size(G,2))];
 Q=[Q zeros(size(Q,1),1)];
 M=[M  zeros(size(M,1),1)];
 M=[M ; zeros(1,size(M,2))];
 M(end,end)=1;
 M(1:ns,end)=const_coef;
end;


end