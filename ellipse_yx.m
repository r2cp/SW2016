function [x y] = ellipse_yx(m,V,c);

% Find X and Y such that
% [(X-mx) (Y-my)]'V[(X-mx) (Y-my)] = c
%
%  Input m: [mx my]';
%        V: a 2 x 2 matrix
%        c: constant
%


%Step 1: Trace out Unit Sphere
n = 500;
x = linspace(-1,1,n)';
y = sqrt(ones(n,1)-x.^2);
x = sqrt(c)*x;
y = sqrt(c)*y;
z1 = [x y]';
z2 = [x -y]';
cv = chol(V);
cvi = inv(cv);
q1 = (cvi*z1 + repmat(m,1,n))';
q2 = (cvi*z2 + repmat(m,1,n))';
x = [q1(:,1) q2(:,1)];
y = [q1(:,2) q2(:,2)];


% Example of Code to plot the ellipse
% plot(x(:,1),y(:,1),'-k','LineWidth',3);
% hold on;
%  plot(x(:,2),y(:,2),'-k','LineWidth',3);
% hold off;

end

