function [ HP_WGHT ] = hp_weights(n,alpha )
% This function computes the HP Weighting matrix for a vector of size n
% using a weighting function alpha
% y_hp = HP_WGHT*y
% This function is a slightly modified version of the program written by
% Adriana Cornea-Madeira 
% and is described below Thm 2.1 in her paper 'The Explicit Formula for the
% Hodrick-Prescott Filter in Finite Sample" (Review of Economics and
% Statistics, Forthcoming 2016).

i = linspace(1,n,n)';
j = linspace(1,n,n)';
T = sin(i*j'*pi/(n+1));
la = 1+alpha*(2-2*cos(i*pi/(n+1))).^2;
La = diag(la);
Lai = diag(1./la);
T = sqrt(2)*T/sqrt(n+1); %matrix of eigenvectors
S = sqrt((n+1)/2);

% computation of matrix K1 with typical element (19)
i = 1:2:n;
j = 1:2:n;
K1 = zeros(n,n);
for k = 1:length(i)
    for m = 1:length(j)
        K1(i(k),j(m)) = (2*sin(i(k)*pi/(n+1))/S-sin(2*i(k)*pi/(n+1))/S)*(2*sin(j(m)*pi/(n+1))/S-sin(2*j(m)*pi/(n+1))/S)/((1+4*alpha*(1-cos(i(k)*pi/(n+1)))^2)*(1+4*alpha*(1-cos(j(m)*pi/(n+1)))^2));
    end
end
h = sum(((2*sin(j*pi/(n+1))/S - sin(2*j*pi/(n+1))/S).^2)./(1+alpha*(2-2*cos(j*pi/(n+1))).^2));
k1 = 2*alpha./(1-2*alpha*h); % the constant k1

% computation of matrix K2 with typical element (19)
i = 2:2:n;
j = 2:2:n;
K2 = zeros(n,n);
for k = 1:length(i)
    for m = 1:length(j)
        K2(i(k),j(m)) = (2*sin(i(k)*pi/(n+1))/S-sin(2*i(k)*pi/(n+1))/S)*(2*sin(j(m)*pi/(n+1))/S-sin(2*j(m)*pi/(n+1))/S)/((1+4*alpha*(1-cos(i(k)*pi/(n+1)))^2)*(1+4*alpha*(1-cos(j(m)*pi/(n+1)))^2));
    end
end
d = sum(((2*sin(j*pi/(n+1))/S - sin(2*j*pi/(n+1))/S).^2)./(1+alpha*(2-2*cos(j*pi/(n+1))).^2));
k2 = 2*alpha./(1-2*alpha*d); % the constant k2

% weights from Theorem 2.1
Finv_new = T*Lai*T + k1*T*K1*T + k2*T*K2*T;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparison with the weights from the HP filter based on the direct
% inversion of the matrix
% p = n:-1:1;
% P = eye(n);
% P = P(p,:);
% Q = diag(2*ones(n,1))+diag(-1*ones(n-1,1),1)+diag(-1*ones(n-1,1),-1);
% g = [-2;1;zeros(n-2,1)];
% A = eye(n)+alpha*Q*Q;
% F = A-alpha*g*g'-alpha*P*g*g'*P;
% Finv_old = inv(F); % Finv_old should be identical to Finv_new

HP_WGHT = Finv_new;

end

