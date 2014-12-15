function reconCost(trueValue, order, ad_matrix, nad_matrix, lost_matrix)
%reconCost
%trueValue, a vector which contain real reconciliation cost
%  MATRIX: a matrix data with each column as Tree dlc for each method

%matrix doit être construit de cette facon: RAxML.ad, TreeFix.ad,
%PolySolver.95.ad, Polysolver*

perf_dup= compute_performance(trueValue.TruePhylo_dup, ad_matrix);
perf_lost= compute_performance(trueValue.TruePhylo_lost, lost_matrix);
perf_nad= compute_performance(trueValue.TruePhylo_nad, nad_matrix);

figure,

ax1=subplot(3,1,1);
suplot(perf_dup, ax1,'app. dup accuracy, precision and sensitivity for each program', order);

ax2=subplot(3,1,2);
suplot(perf_lost, ax2,'lost accuracy, precision and sensitivity for each program', order);

ax3=subplot(3,1,3);
suplot(perf_nad, ax3,'non app. dup accuracy, precision and sensitivity for each program', order);

end

function suplot(perf, ax, titre, order)

h=bar(ax,perf);

ybuff=5;
for i=1:length(h)
    XDATA=get(get(h(i),'Children'),'XData');
    YDATA=get(get(h(i),'Children'),'YData');
    for j=1:size(XDATA,2)
        x=XDATA(1,j)+(XDATA(3,j)-XDATA(1,j))/2;
        y=YDATA(2,j)+ybuff;
        t=num2str(YDATA(2,j),3);
        if(strcmp(t, '0'))
             text(x,y+ybuff*5,'*','Color','r','HorizontalAlignment','left', 'FontSize', 15, 'fontWeight','bold')
        else
            text(x,y,t,'Color','k','HorizontalAlignment','left','Rotation',90, 'FontSize', 10)
        end
    end
end

ylabel('%');
ylim([0, 145])
xlim([0.5, numel(order)+2])
title(titre);
legend({'Accuracy','Precision', 'Sensibility'});
set(ax,'XTickLabel',order,'TickLength',[0.001 0])


end

function performance= compute_performance(true_value, matrix)

n_element=size(matrix, 2);
true_positive=zeros(1, n_element); % duplication exist and program find the exact number
true_negative=true_positive; % duplication is 0 and program find 0 duplication
sensi_deno=true_positive;% duplication inferred < true duplication 
preci_deno=true_positive; % duplication inferred > true duplication

%false_positive can be changed to duplication is 0 and program infered
%suplication
%so is false negative which can be changed to : duplication exist and
%program infered 0 duplication
for i=1:n_element
    true_positive(i)=nnz((true_value==matrix(:,i) & true_value>0));
    true_negative(i)=nnz((true_value==matrix(:,i) & true_value==0));
   
    sensi_deno(i)=nnz(true_value>0);
    preci_deno(i)=nnz(matrix(:,i)>0);
end


acc= (true_positive + true_negative)*100/size(matrix,1);

sen= (true_positive*100)./sensi_deno;

pre= (true_positive*100)./preci_deno;

performance = [acc', pre', sen'];
end

