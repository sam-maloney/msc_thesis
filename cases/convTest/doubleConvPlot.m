eps = '0.1';

% N = [25 50 100 200 400 800 1600 3200 6400];
N = [10 30 90 270 810 2430 7290];

ref_rho = load(['firstOrderfalse' eps '/ref_rho'], '-ascii');
ref_rho = sortrows(ref_rho);

limiter = {'firstOrder','Venkatakrishnan'};
order = {'First Order','Second Order'};

for k = 1:2

    err_unb = zeros(size(N));
    err_bal = zeros(size(N));

    for i = 1:length(N)
        rho = load([limiter{k} 'false' eps '/' num2str(N(i)) '_rho'], '-ascii');
        interp_ref = interp1(ref_rho(:,1), ref_rho(:,2), rho(:,1), 'nearest');
        err_unb(i) = norm(rho(:,2) - interp_ref, 1)/N(i);

        rho = load([limiter{k} 'true' eps '/' num2str(N(i)) '_rho'], '-ascii');
    %     interp_ref = interp1(ref_rho(:,1), ref_rho(:,2), rho(:,1), 'spline');
        err_bal(i) = norm(rho(:,2) - interp_ref, 1)/N(i);
    end

    dx = 4./N;
    
    subplot(1,2,k);

    %%%%% plot error norm magnitudes %%%%%
    yyaxis left;
    mSize = 20;
    loglog(dx,err_unb,'--d');
    hold on;
    loglog(dx,err_bal,'--.','MarkerSize',mSize);
    hold off;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    x_lims = [min(dx)/sqrt(2) max(dx)*sqrt(2)];
    xlim(x_lims);
    xticks([10^-3 10^-2 10^-1]);
    ylim([0.4*min([err_unb err_bal]) 2*max([err_unb err_bal])]);
%     ylim([0.25*min([err_unb err_bal]) 4*max([err_unb err_bal])]); % eps=0,0.0001
    yticks([10^-10 10^-8 10^-6 10^-4 10^-2]); % eps=0.0001
    yticks([10^-14 10^-10 10^-6 10^-2]); % eps=0
    set(gca,'XDir','reverse');
    xlabel('\Deltax');
    ylabel('Normalized L1 Error Norm');
    title(order{k});

    %%%%% 2-axis order plot %%%%%
    hold on;
    yyaxis right;

    l_dx = log(dx);
    dx2 = dx(2:end)*sqrt(2);

    l_err = log(err_unb);
    order_unb = (l_err(2:end)-l_err(1:end-1))./(l_dx(2:end)-l_dx(1:end-1));
    semilogx(dx2,order_unb,'d');

    l_err = log(err_bal);
    order_bal = (l_err(2:end)-l_err(1:end-1))./(l_dx(2:end)-l_dx(1:end-1));
    semilogx(dx2,order_bal,'.','MarkerSize',mSize);

    ylim([-1 3]);
    yticks([0 1 2]);
    ylabel('Order of Accuracy');
    plot(x_lims,[0 0],':',x_lims,[1 1],':',x_lims,[2 2],':');

    yyaxis left;
    
%     if ( k == 1 )
        h = zeros(1,2);
        h(1) = plot(0,0,'dk');
        h(2) = plot(0,0,'.k','MarkerSize',mSize);
        legend(h,'Unbalanced','Well-Balanced','location','SouthWest');
%     end
    
    hold off;

%     legend('Unbalanced Err_1','Balanced Err_1','Unbalanced order','Balanced order');
    
end