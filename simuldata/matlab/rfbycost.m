function rfbycost (rf_value, maxrf_value, true_recon, order, roundValue, lim)

if nargin <5, roundValue=1; lim=inf; end

nnz(true_recon>lim)
numel(true_recon)

reconcost=round(true_recon/roundValue)*roundValue;
recon_uniq= unique(reconcost);
recon_uniq=recon_uniq(recon_uniq<lim);

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
%ax2=subplot(2,1,1);
ax2= gca;
plot(recon_uniq, rf_accu_plot)
xlim([0, max(recon_uniq)+1.55]);
xlabel('Reconciliation cost');
ylabel('Topology Accuracy (RF=0)');
legend2=legend(order);
set(legend2,'FontSize',fsize);
axt1=title('RAxML, ProfileNJ, TreeFix topology accuracy for increasing number of evolutionary event');

set(ax2,'box', 'off', 'FontSize', fsize, 'LineWidth', 1.2,'FontName', 'Helvetica', 'TickDir', 'out')
set(axt1, 'FontWeight', 'bold','FontSize',fsize+1, 'FontName', 'Helvetica');

% ax3=subplot(2,1,2);
% area(recon_uniq, rf_mean_plot)
% xlim([0, max(recon_uniq)+1.55])
% xlabel('Reconciliation cost');
% ylabel('mean RF');
% colormap jet;
% l2=legend(order);
% set(l2,'FontSize',fsize);
% axt2=title('mean RF distance between true tree and RAxML, ProfileNJ, TreeFix for increasing number of evolutionary event');
% 
% set(ax3,'box', 'off', 'FontSize', fsize,'LineWidth', 1.2, 'FontName', 'Helvetica', 'TickDir', 'out')
% set(axt2, 'FontWeight', 'bold','FontSize',fsize+1, 'FontName', 'Helvetica');

set_figures(h, fsize)

% set(gcf,'NextPlot','add');
% axes;
% 
% set(gca,'Visible','off');
% set(axt1,'Visible','on');
% set(axt2,'Visible','on');

end

function set_figures(fig, fsize)
    set(fig,'units','normalized','outerposition',[0 0 1 0.95])
    set(findall(gcf,'type','text'),'FontSize',fsize, 'FontName', 'AvantGarde')
    set(fig, 'PaperPositionMode', 'auto');
    set(fig,'InvertHardcopy','on');
    set(fig,'PaperUnits', 'inches');
end
