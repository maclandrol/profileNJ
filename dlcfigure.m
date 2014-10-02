function dlcfigure(matrix, order)
%RFfigure(MATRIX)
%  MATRIX: a matrix data with each column as Tree dlc for each method

figure1 = figure;
shape={'o', '+', 's', 'p', 'v', '*', 'x', 'd'};
colors= distinguishable_colors(numel(order));
for i=1:size(matrix, 2)
   scatter(1:size(matrix, 1), matrix(:, i),shape{i}, 'MarkerEdgeColor',colors(i, :) , 'MarkerFaceColor', colors(i, :),'LineWidth',0.5);
   hold on;
end
yt=cellstr(get(gca, 'YtickLabel'));
yt(strcmp('1', yt))={'Correct cost -->'};
set(gca, 'YtickLabel', yt);
ylabel('computed dlc/ expected dlc');
xlabel('tree');
legend(order);

title('Comparison of duplication-lost cost between different contraction threshold');

plotdlc(matrix, order)

end


function plotdlc(matrix, order)

%RFfigure(MATRIX)
%  MATRIX: a matrix data with each column as Tree dlc for each method
figure2 = figure;
colors= distinguishable_colors(numel(order));
for i=1:size(matrix, 2)
    plot(matrix(:, i), 'color', colors(i,:));
    hold on
end
legend(order)

end