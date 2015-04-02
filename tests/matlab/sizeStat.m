function sizeStat(treesize, njtree, ad, lost, nad, recon, time, rf, stat_order, time_label,datatype, roundValue)
%Compute statistique based on tree size
% true_dlc est une matrice 

if nargin ==11, roundValue=1; end
treesize=round(treesize/roundValue)*roundValue;
unique_size= unique(treesize);
unique_size=unique_size(unique_size<=60);
dup_size=zeros(numel(unique_size), numel(stat_order));
recon_size= zeros(numel(unique_size), numel(stat_order));
lost_size=zeros(numel(unique_size), numel(stat_order));
rf_size=zeros(numel(unique_size), numel(stat_order));
time_size=zeros(numel(unique_size), numel(time_label));

dup=ad+nad;
dup_accuracy = accuracy(dup, njtree.TruePhylo_dup+njtree.TruePhylo_nad);
lost_accuracy= accuracy(lost, njtree.TruePhylo_lost);
recon_accuracy = accuracy(recon, njtree.TruePhylo_dlc);

for i=1:numel(unique_size)
    ind_i=(treesize==unique_size(i));
    dup_size(i, :)= (sum(dup_accuracy(ind_i, :), 1)*100.0)./nnz(ind_i);
    lost_size(i, :)= (sum(lost_accuracy(ind_i, :), 1)*100.0)./nnz(ind_i);
    recon_size(i, :)= (sum(recon_accuracy(ind_i, :), 1)*100.0)./nnz(ind_i);
    rf_size(i, :) = sum(rf(ind_i, :)==0, 1)*100.0./nnz(ind_i);
    time_size(i, :)= mean(time(ind_i, :), 1);
end

fsize=16;
h1=figure;
plot(unique_size, dup_size);
ylim([0, 130])
legend(stat_order,'FontSize', fsize);
xlabel('number of genes per tree', 'FontSize', fsize);
title(['Accuracy of inferred gene duplication for gene trees of increasing size in simulated ', datatype, ' dataset'],'FontWeight', 'bold','FontSize', fsize+1);
ylabel('duplication accuracy (%)', 'FontSize', fsize);
set_figures(h1, gca)


h2=figure; 
plot(unique_size, lost_size);
ylim([0, 130])
legend(stat_order,'FontSize', fsize);
xlabel('number of genes per tree', 'FontSize', fsize);
title(['Accuracy of inferred gene loss for gene trees of increasing size in simulated ', datatype,' dataset'],'FontWeight', 'bold','FontSize', fsize+1);
ylabel('loss accuracy (%)', 'FontSize', fsize);
set_figures(h2, gca)

h3=figure; 
plot(unique_size, recon_size);
ylim([0, 130])
legend(stat_order,'FontSize', fsize);
xlabel('number of genes per tree', 'FontSize', fsize);
title(['Accuracy of reconciliation cost for gene trees of increasing size in simulated ', datatype,' dataset'],'FontWeight', 'bold','FontSize', fsize+1);
ylabel('reconciliation cost accuracy (%)', 'FontSize', fsize);
set_figures(h3, gca)

h4=figure; 
plot(unique_size, rf_size);
ylim([0, 130])
legend(stat_order, 'FontSize', fsize);
xlabel('number of genes per tree', 'FontSize', fsize);
ylabel('Topologies correct (%)', 'FontSize', fsize);
title(['Accuracy of the inferred topology for gene trees of increasing size in simulated ', datatype,' dataset'],'FontWeight', 'bold','FontSize', fsize+1);
set_figures(h4, gca)

h5=figure; 
plot(unique_size, time_size);
ylim([0, max(max(time_size))*1.5])
legend(time_label, 'FontSize', fsize-1);
xlabel('number of genes per tree', 'FontSize', fsize);
ylabel('time (s)', 'FontSize', fsize);
title(['Runtime for gene trees of increasing size in simulated ',datatype, ' dataset'],'FontWeight', 'bold','FontSize', fsize+1);
set_figures(h5, gca)

end


function acc=accuracy(matrix, trueval)
    acc=bsxfun(@eq,matrix,trueval);
end


function set_figures(fig, axis)
    set(fig,'units','normalized','outerposition',[0 0 1 0.95])
    set(fig, 'PaperPositionMode', 'auto');
    set(fig,'InvertHardcopy','on');
    set(fig,'PaperUnits', 'inches');
    set(axis, 'TickDir', 'out',   'TickLength', [.01 .0001], 'box', 'off');
end


