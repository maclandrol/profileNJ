function timesol(times, npsol, labels, roundValue)

if nargin ==3, roundValue=1; end
npsol=round(npsol/roundValue)*roundValue;
uniq_sol= unique(npsol);

% Create figure
mean_times= zeros(size(uniq_sol, 1), numel(labels));

for i=1:numel(uniq_sol)
    ind=(npsol==uniq_sol(i));
    mean_times(i,:)= mean(times(ind, :),1);
    times(ind, :)

end

h=figure;
fsize=14;
plot(uniq_sol, mean_times)
xlabel('number of solutions');
ylabel('times (s)');
legend(labels, 'location', 'northwest')
mean_times(6,:)
set(gca, 'box','off', 'TickDir', 'out');

set_figures(h, fsize)

set(gcf,'NextPlot','add');
axes;
ht = title('profileNJ runtime as a function of the number of solutions', 'FontWeight', 'bold','FontSize',fsize+1, 'FontName', 'Helvetica');
set(gca,'Visible', 'off');
set(ht,'Visible','on');

end

function set_figures(fig, fsize)
    set(fig,'units','normalized','outerposition',[0 0 1 0.95])
    set(findall(gcf,'type','text'),'FontSize',fsize, 'FontName', 'AvantGarde')
    set(fig, 'PaperPositionMode', 'auto');
    set(fig,'InvertHardcopy','on');
    set(fig,'PaperUnits', 'inches');
end
