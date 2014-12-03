function RFfigure(matrix, order)
%RFfigure(MATRIX)
%  MATRIX: a matrix data with each column as Tree RF for each method


% Create figure
figure1 = figure;

set(gcf,'units','normalized','outerposition',[0 0 1 0.95])
fsize= 16;
set(findall(figure1,'type','text'),'FontSize',fsize, 'FontName', 'AvantGarde')
set(figure1, 'PaperPositionMode', 'auto');
set(figure1,'InvertHardcopy','on');
set(figure1,'PaperUnits', 'inches');

% Create axes
binrange=0:2:max(matrix(:));
bincount= histc(matrix, binrange);
bincount= bincount*100.0/size(matrix,1);

% Create multiple lines using matrix input to bar
bar(binrange, bincount,'grouped');
set(gca, 'YTick',[0 25 50 75 100],...
    'TickLength',[0.002 0],...
    'Xtick', binrange,...
    'TickDir', 'out', ...
    'FontSize',fsize,...
    'LineWidth', 1.2,...
    'FontName','Arial', 'box','off');
xlim([-2,29]);
% Create xlabel
xlabel('RF distance')

% Create ylabel
ylabel('% of trees')

% Create legend
legend1 = legend(order);
set(legend1,'FontSize',18);

title('Topology accuracy of RAxMl, ProfileNJ and TreeFix trees mesured by the RF distance with the simulated tree', 'FontSize', 18, 'FontWeight' , 'bold');
%