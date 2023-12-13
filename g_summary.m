function [g_prior_post] = g_summary(g_prior,g_draws)
% Discrete prior for g is in g_prior
%  first column = values, second column = prior probs
%  g_draws is (n_draws) x (n_s) with draws of g from posterior for n_s
%  series
%  On out g_prior_post has 
%   Values of g in column 1
%   Prior in column 2
%   Posterior in columns 3-n_s+1

 n_g = size(g_prior,1);
 n_s = size(g_draws,2);
 g_prior_post = [g_prior ones(n_g,n_s)];
 for i = 1:n_g;
     for j = 1:n_s;
         tmp = g_draws(:,j);
         g_prior_post(i,2+j) = mean(tmp == g_prior(i,1));
     end;
 end;

end

