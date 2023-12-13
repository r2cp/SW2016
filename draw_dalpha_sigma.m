function [sigma_draw] = draw_dalpha_sigma(dalpha,nu_prior_alpha, s2_prior_alpha)

T = size(dalpha,1);
n = size(dalpha,2);
SSR_mat = sum(dalpha.^2);
SSR_prior = nu_prior_alpha*s2_prior_alpha;
nu = T + nu_prior_alpha;
a = nu/2;
sigma_draw = NaN*zeros(n,1);
for i = 1:n;
  ssr = SSR_mat(1,i) + SSR_prior;
  b = 2/ssr;
  var_dalpha = 1./gamrnd(a,b);
  sigma_draw(i) = sqrt(var_dalpha);
end;

 
end