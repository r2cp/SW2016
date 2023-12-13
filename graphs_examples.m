% Example of Graphics Programs
% For future reference
clear all;
figdir = '/Users/mwatson/Dropbox/';
%
% Time Series Plots
[dnobs,calvec,calds] = calendar_make(1950,2013,1,4,4);
x = cumsum(randn(dnobs,1));
xstd = std(x);
x_upper = x + xstd;
x_lower = x - xstd;
% -------------------- Line graphs ----------------- 
% Some parameters 
% 'LineWidth',xx     where xx is a number
% Line Style
%    '-' solid line
%    '--' dashed lind
%    ':' dotted line
%    '-.' dash dot line
%    'none'
% Colors
%     y yellow
%     m magenta
%     c cyan
%     r red
%     g green
%     b blue
%     w white
%     k black
% FontWeight 
% 'normal'
% 'bold'
% 'light'
% 'demi'

% Plot first to get axis
plot(calvec,[x_lower,x_upper]);
v = axis;  % This puts the axis dimensions in v
axis([1960 2010 v(3) v(4)]);  % Set axis (v(3) and v(4) are the verticl axis)
str_temp = 'Fig';
save_g_str = [figdir str_temp '_1.jpg'];
xlabel('Date','FontSize',12,'FontWeight','normal','Color','k'); 
ylabel('Y-Value','FontSize',12,'FontWeight','normal','Color','k');
  plot(calvec,x,': k','LineWidth',2);
hold on;  
  plot(calvec,x_upper,'-- b','LineWidth',1);
  plot(calvec,x_lower,'- r','LineWidth',1);
hold off;
saveas(gcf,save_g_str);
%
% Scatter plot
%
x = randn(30,1);
y = x + randn(30,1);
plot(x,y,'ok');  % black circles
plot(x,y,'ok','MarkerSize',10); % larger 
plot(x,y,'ok','MarkerSize',10,'MarkerFaceColor','r'); % Fill with red

           