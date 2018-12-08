% N = [25 50 100 200 400 800 1600 3200 6400];
N = [10 30 90 270 810 2430 7290];
err_1 = zeros(size(N));
err_Inf = zeros(size(N));
max_diff_val = zeros(length(N),1);
max_diff_idx = zeros(length(N),1);

ref_rho = load('ref_rho', '-ascii');
ref_rho = sortrows(ref_rho);

for i = 1:length(N)
    rho = load([num2str(N(i)) '_rho'], '-ascii');
    interp_ref = interp1(ref_rho(:,1), ref_rho(:,2), rho(:,1), 'spline');
    diff = rho(:,2) - interp_ref;
    err_1(i) = norm(diff, 1)/N(i);
    err_Inf(i) = norm(diff, Inf);
    [max_diff_val(i), max_diff_idx(i)] = max(abs(diff));
end

max_diff = [max_diff_val, max_diff_idx];

dx = 4./N;

max_diff(:,3) = dx'.*(max_diff(:,2) - 0.5) - 2;

%%%%% plot error norm magnitudes %%%%%
yyaxis left;
mSize = 17;
loglog(dx,err_1,'.','MarkerSize',mSize);
hold on;
loglog(dx,err_Inf,'d');
hold off;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% linear fit order plot %%%%%
% hold on;
% 
% % linear fit for L-1 norm
% % fit_range = 2:length(dx)-5;
% fit_range = 8:length(dx);
% B = polyfit(log(dx(fit_range)) ,log(err_1(fit_range)), 1);
% y = exp(polyval(B, log([dx(1) dx(end)])));
% 
% loglog([dx(1) dx(end)],y,'--b');
% text(dx(end-3)*0.9,err_1(end-3),['     O(\Deltax^{' num2str(B(1)) '})']);
% 
% % linear fit for L-Inf norm
% % fit_range = 8:length(dx);
% fit_range = 2:length(dx)-5;
% B = polyfit(log(dx(fit_range)) ,log(err_Inf(fit_range)), 1);
% y = exp(polyval(B, log([dx(1) dx(end)])));
% 
% loglog([dx(1) dx(end)],y,'--r');
% text(dx(end-3)*0.9,err_Inf(end-3),['     O(\Deltax^{' num2str(B(1)) '})']);
% 
% hold off;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x_lims = [min(dx)/sqrt(2) max(dx)*sqrt(2)];
xlim(x_lims);
ylim([0.5*min([err_1 err_Inf]) 2*max([err_1 err_Inf])]);
% yticks([10^-6 10^-5 10^-4]);
set(gca,'XDir','reverse');
xlabel('\Deltax');
ylabel('Magnitude of error norm');
title('Convergence Study');

%%%%% 2-axis order plot %%%%%
hold on;
yyaxis right;

l_dx = log(dx);
dx2 = dx(2:end)*sqrt(2);

l_err = log(err_1);
order_1 = (l_err(2:end)-l_err(1:end-1))./(l_dx(2:end)-l_dx(1:end-1));
semilogx(dx2,order_1,'--.','MarkerSize',mSize);

l_err = log(err_Inf);
order_Inf = (l_err(2:end)-l_err(1:end-1))./(l_dx(2:end)-l_dx(1:end-1));
semilogx(dx2,order_Inf,'--d');

ylim([-1 3]);
yticks([-1 0 1 2 3]);
ylabel('Order of convergence');
plot(x_lims,[0 0],':',x_lims,[1 1],':',x_lims,[2 2],':');

yyaxis left;
hold off;

% legend('L_1 magnitude','L_{Inf} magnitude','L_1 order','L_{Inf} order');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%