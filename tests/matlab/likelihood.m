function likelihood(lkl_matrix, p_value_matrix, rf_matrix, max_rf_matrix, order)

rf_ratio= rf_matrix./max_rf_matrix;
colors= distinguishable_colors(numel(order));
h1=figure;
h2=figure;
fsize=15;
set_figures(h1,fsize)
set_figures(h2,fsize)

shapes={'+', 'o', '^', 's', 'd', '*'};
cols={'r', 'b', 'g', 'k', 'y', 'c'};

for i=1:numel(order)-1
    %Figure1 setting
    figure(h1)
    ax1=subplot(1,3,i);
    semilogy(rf_ratio(:,i), lkl_matrix(:,i), 'x', 'color', colors(i,:))
    %scatter( 'x', 'MarkerEdgeColor', );
    ylabel('Log Likelihood', 'FontSize', 14);
    xlabel('RF distance', 'FontSize', 14);
    set(ax1, 'FontName'    , 'Helvetica',...
    'Box'         , 'off'     , ...
    'TickDir'     , 'out'     , ...
    'TickLength'  , [.01 .001] , ...
    'LineWidth'   , 1.5        , ...
    'FontSize'    ,fsize);
    xlim(ax1, [-0.09, 1]);
    title(sprintf(['Likelihood of ', order{i}, ' trees according to \nthe RF distance with the true tree']));
    
    %Figure 2 setting
    figure(h2)
    ax2=subplot(1,3,i);
    scatter(rf_ratio(:,i), p_value_matrix(:,i), 'x', 'MarkerEdgeColor',colors(i,:) );
    xlim(ax2, [-0.09, 1]);
    ylabel('SH-test pval', 'FontSize', 14);
    xlabel('RF distance', 'FontSize', 14);
    set(ax2, ...
    'FontName'    , 'Helvetica',...
    'Box'         , 'off'     , ...
    'TickDir'     , 'out'     , ...
    'TickLength'  , [.01 .001] , ...
    'YMinorTick'  , 'on'      , ...
    'LineWidth'   , 1.5        , ...
    'FontSize'    ,fsize);
    title(sprintf(['SH-test p-val of ', order{i}, ' trees according to \nthe RF distance with the true tree']));
    
end



%Stack likelihood on the same plot
h3=figure;
set_figures(h3,fsize)

for i=1:numel(order)-1
    semilogy(rf_ratio(:,i), lkl_matrix(:,i), [shapes{i} cols{i}])
    hold on
end
 
set(gca, ...
    'FontName'    , 'Helvetica',...
    'Box'         , 'off'     , ...
    'TickDir'     , 'out'     , ...
    'TickLength'  , [.01 .001] , ...
    'YMinorTick'  , 'on'      , ...
    'LineWidth'   , 1.5        , ...
    'FontSize'    ,fsize);
xlim([-0.09, 1]);
ylabel('Log-Likelihood', 'FontSize', 14);
xlabel('RF distance', 'FontSize', 14);
legend(order{1:end-1}, 'Location', 'NorthEast');
title(sprintf('Tree likelihood according to \nthe RF distance with the true tree for each program'));
    

% Box plot setting
h4=figure;

boxplot(p_value_matrix, 'labels', order, 'notch','on', 'colors', distinguishable_colors(numel(order)), 'boxstyle', 'outline')
ylabel('AU test p-value');
title(sprintf('Comparision of AU test p-value between the true tree \nand tree returned by RAxML, profileNJ and TreeFix'), 'FontName', 'Helvetica','FontSize', fsize, 'FontWeight', 'bold');
set(gca, 'YTick', [0, 0.05, 0.2, 0.4, 0.6, 0.8, 1]);
set(gca, 'YTickLabel', {'0','a=0.05', '0.2', '0.4', '0.6', '0.8', '1'});
hold on
plotp=plot(0:1000, 0.05*ones(1001,1), '--');
set(plotp, 'color', [0.5,0.5,0.5]);
set(findall(h4,'type','text'),'FontSize',14);
set(h4, 'PaperPositionMode', 'auto');
set(h4,'InvertHardcopy','on');
set(gca,'TickDir', 'out', 'FontSize', 14, 'Box', 'off','TickLength'  , [.01 .01]);

set(gcf,'units','normalized','outerposition',[0 0 0.7 0.6])
set(gcf, 'PaperPositionMode', 'auto');

end

function set_figures(fig, fsize)
    set(fig,'units','normalized','outerposition',[0 0 1 0.95])
    set(findall(gcf,'type','text'),'FontSize',fsize, 'FontName', 'AvantGarde')
    set(fig, 'PaperPositionMode', 'auto');
    set(fig,'InvertHardcopy','on');
    set(fig,'PaperUnits', 'inches');
end