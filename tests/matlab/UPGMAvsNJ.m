function UPGMAvsNJ(UPGMA, NJ, order)

figure, 

upgma_accuracy= mean(UPGMA==0, 1)*100;
nj_accuracy= mean(NJ==0, 1)*100;

rf_accuracy=[upgma_accuracy', nj_accuracy'];

h=bar(rf_accuracy, 0.9)
ybuff=5;
for i=1:length(h)
    XDATA=get(get(h(i),'Children'),'XData');
    YDATA=get(get(h(i),'Children'),'YData');
    for j=1:size(XDATA,2)
        x=XDATA(1,j)+(XDATA(3,j)-XDATA(1,j))/2;
        y=YDATA(2,j)+ybuff;
        t=num2str(YDATA(2,j),3);
        text(x,y,t,'Color','k','HorizontalAlignment','left','Rotation',90, 'FontSize', 10)
    end
end

set(h(1),  'FaceColor',[61 219 172]./255)
set(h(2),  'FaceColor', [208, 89, 89]./255)%'g','EdgeColor','k','LineWidth',3)

ylabel('% of correct topologies');
xlabel('Contraction threshold \lambda');
set(gca,'box','off')

ylim([0, 130]);
xlim([0.5, numel(order)+0.8]);
set(gca, 'XTickLabel', order);
title('Comparision of UPGMA and NJ clustering method in PolytomySolver');
legend({'UPGMA', 'NJ'});

