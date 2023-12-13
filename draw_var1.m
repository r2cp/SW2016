proc(2) = draw_var1(x,r_pt,gam,var_min,r_sig2,r_m1,r_m2,r_p);

@ Construct draw of variance @
local vague, n,lnres2,tmp,ir,mut,adraw1,edraw;
local f1,f2,fe,vardraw,r_sig;
local ye,p1t,p2t,x1t,x2t,x3_draw,x1,x2,x3,p1,p2,p3,p2i,k,t,h,utmp,x3mean,chol_p;
r_sig=sqrt(r_sig2);

 vague = 1.0e+6;  @ large number @    
 n=rows(x);
 @ Construct estimate of var_e @
 lnres2=ln(x.^2);
     
 @ -- Step 1 -- initial draws of Indicator Variables -- @
 tmp=rndu(n,1);
 ir=tmp .<= r_pt;

 @ -- Step 2; compute system parameters given indicators -- @
 mut = (ir*r_m1) + ((1-ir)*r_m2);

 @ Kalman Filtering @
 ye=lnres2-mut;
 p1t=zeros(n+1,1);
 p2t=zeros(n+1,1);
 x1t=zeros(n+1,1);
 x2t=zeros(n+1,1);
 x3_draw=zeros(n+1,1);

 @ -- Compute Covariance Matrix  -- @
 x1=0;
 p1=vague;
 x1t[1,.]=x1; 
 p1t[1,.]=p1;
 for t (1,n,1);
 	x2=x1;
 	p2=p1+gam;
 	h=p2+r_sig2;
 	k=p2/h;
 	p1=p2-k*p2;
 	x1=x2+k*(ye[t]-x2);
 	p1t[t+1]=p1;
 	p2t[t+1]=p2;
 	x1t[t+1]=x1;
 	x2t[t+1]=x2;
 endfor;
 utmp=rndn(n,1);
 x3mean=x1;
 p3=p1;
 chol_p=sqrt(p3);
 x3=x3mean+chol_p*utmp[n];
 x3_draw[n+1,.]=x3';
 for t (n,2,-1);
  x2=x2t[t+1];
  p2=p2t[t+1];
  x1=x1t[t];
  p1=p1t[t];
  p2i=1/p2;
  k=p1*p2i;
  x3mean=x1+k*(x3-x2);
  p3=p1-k*p1;
  chol_p=sqrt(p3);
  x3=x3mean+chol_p*utmp[t];
  x3_draw[t,.]=x3;
 endfor;
 adraw1=x3_draw[2:n+1];
 edraw=lnres2-adraw1;

 @ -- Compute Mixture Probabilities -- @
  @ -- e shocks -- (Note sigma is the same) -- @
  f1=exp(   (-0.5)* (((edraw-r_m1)./r_sig).^2)  );
  f2=exp(   (-0.5)* (((edraw-r_m2)./r_sig).^2)  );
  fe= r_p*f1 + (1-r_p)*f2;
  r_pt=(r_p*f1)./fe;
  
   @ -- Compute Variance of AR Shocks and Variance of Annual Difs -- @
  vardraw = exp(adraw1); @ Variance Draw @;
  
 @ -- Impose minimum value --- @
 tmp=var_min*ones(rows(vardraw),1);
 vardraw=maxc((vardraw~tmp)'); 
 retp(vardraw,r_pt);
endp;