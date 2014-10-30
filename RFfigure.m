function RFfigure(matrix, order)
%RFfigure(MATRIX)
%  MATRIX: a matrix data with each column as Tree RF for each method


% Create figure
figure1 = figure;

% Create axes
binrange=0:2:max(matrix(:));
bincount= histc(matrix, binrange);
bincount= bincount*100.0/size(matrix,1);

% Create multiple lines using matrix input to bar
bar(binrange, bincount, 'grouped');
set(gca, 'YTick',[0 25 50 75 100],...
    'TickLength',[0.001 0],...
    'Xtick', binrange,...
    'FontSize',14,...
    'FontName','Arial', 'box','off');

% Create xlabel
xlabel('RF distance','FontSize',15);

% Create ylabel
ylabel('% of trees','FontSize',15);

% Create legend
legend1 = legend(order);
set(legend1,'FontSize',18);
title('RF comparision between TreeSolver, TreeFix and RAxML method')
%