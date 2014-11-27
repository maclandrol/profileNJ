function reconcost(trueValue, order, dup_matrix, loss_matrix, recon_matrix)
%reconCost
%trueValue, a vector which contain real reconciliation cost, in this order:
%ad, lost, recon
%  MATRIX: a matrix data with each column as Tree dlc for each method

%matrix doit être construit de cette facon: RAxML.ad, TreeFix.ad,
%PolySolver.95.ad, Polysolver*

n_el= size(dup_matrix, 1);
dupaccuracy= 100*sum(bsxfun(@eq, dup_matrix, trueValue(:,1)))./n_el;
lossaccuracy= 100*sum(bsxfun(@eq, loss_matrix, trueValue(:,2)))./n_el;
reconaccuracy= 100*sum(bsxfun(@eq, recon_matrix, trueValue(:,3)))./n_el;
perf_recon = [dupaccuracy', lossaccuracy', reconaccuracy'];
figure,
set(gcf,'units','normalized','outerposition',[0 0 1 0.95])
fsize= 15;
set(findall(gcf,'type','text'),'FontSize',fsize, 'FontName', 'AvantGarde')
rplot(perf_recon, gca, order, fsize);
set(gcf, 'PaperPositionMode', 'auto');
set(gcf,'InvertHardcopy','on');
set(gcf,'PaperUnits', 'inches');
title('Duplication, Loss and Reconciliation accuracy for each program', 'FontSize', 16, 'FontWeight' , 'bold');

end

function rplot(perf, ax, order, fsize)

h=bar(ax,perf);

ybuff=3;
for i=1:length(h)
    XDATA=get(get(h(i),'Children'),'XData');
    YDATA=get(get(h(i),'Children'),'YData');
    for j=1:size(XDATA,2)
        x=XDATA(1,j)+(XDATA(3,j)-XDATA(1,j))/2;
        y=YDATA(2,j)+ybuff;
        t=num2str(YDATA(2,j),3);
        if(strcmp(t, '0'))
             text(x,y+ybuff*3,'*','Color','r','HorizontalAlignment','left', 'FontSize', fsize, 'FontWeight','bold');
        else
            text(x,y,t,'Color','k','HorizontalAlignment','left','Rotation',90, 'FontSize',fsize-2,'FontName', 'Helvetica' );
        end
    end
end

set(h(1),  'FaceColor',[94 179 86]./255)
set(h(2),  'FaceColor', [210, 44, 44]./255)%'g','EdgeColor','k','LineWidth',
set(h(3),  'FaceColor',[76 90 181]./255)

ylabel('% Accuracy', 'FontName', 'Helvetica', 'FontSize',fsize);
ylim([0, 135])
xlim([0.5, numel(order)+1.5])
legend({'Duplication','Loss', 'Reconciliation'});

set(ax, ...
  'FontName'    , 'Helvetica',...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.01 .001] , ...
  'XTickLabel'  , order     , ...
  'YMinorTick'  , 'on'      , ...
  'LineWidth'   , 1.5        , ...
  'FontSize'    ,fsize);


end

