function rfbycost (rf_value, maxrf_value, true_recon, order, roundValue)

if nargin ==4, roundValue=1; end
reconcost=round(true_recon/roundValue)*roundValue;
recon_uniq= unique(reconcost);

% Create figure
rf_value=rf_value./maxrf_value;
rf_accu_plot= zeros(size(recon_uniq, 1), numel(order));
rf_mean_plot= zeros(size(recon_uniq, 1), numel(order));

for i=1:numel(recon_uniq)
    ind=(reconcost==recon_uniq(i));
    rfs= rf_value(ind, :);
    acc=sum(rfs==0, 1)/size(rfs, 1);
    moy= mean(rfs,1);
    rf_accu_plot(i, :)= acc;
    rf_mean_plot(i, :)= moy;
end

h=figure;
fsize=14;
ax1=subplot(2,1,1);
plot(recon_uniq, rf_accu_plot)
xlim([0, max(recon_uniq)+1.55]);
xlabel('Reconciliation cost');
ylabel('Topology Accuracy (RF=0)');
legend(order)


ax2=subplot(2,1,2);
area(recon_uniq, rf_mean_plot)
xlim([0, max(recon_uniq)+1.55])
xlabel('Reconciliation cost');
ylabel('mean RF');
colormap jet;
legend(order)
set(ax1,'box', 'off', 'FontSize', fsize, 'FontName', 'Helvetica', 'TickDir', 'out')
set(ax2,'box', 'off', 'FontSize', fsize, 'FontName', 'Helvetica', 'TickDir', 'out')

set_figures(h, fsize)

set(gcf,'NextPlot','add');
axes;
ht = title('RAxML, ProfileNJ, TreeFix accuracy for increasing number of evolutionary event', 'FontWeight', 'bold','FontSize',fsize+1, 'FontName', 'Helvetica');
set(gca,'Visible','off');
set(ht,'Visible','on');

end

function set_figures(fig, fsize)
    set(fig,'units','normalized','outerposition',[0 0 1 0.95])
    set(findall(gcf,'type','text'),'FontSize',fsize, 'FontName', 'AvantGarde')
    set(fig, 'PaperPositionMode', 'auto');
    set(fig,'InvertHardcopy','on');
    set(fig,'PaperUnits', 'inches');
end
