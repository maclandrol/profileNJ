function sizeStat(treesize, njtree, dup, lost, nad, time, rf, stat_order, roundValue)
%Compute statistique based on tree size
% true_dlc est une matrice 

if nargin ==8, roundValue=1; end
treesize=round(treesize/roundValue)*roundValue;
unique_size= unique(treesize);
dup_size=zeros(numel(unique_size), numel(stat_order));
lost_size=zeros(numel(unique_size), numel(stat_order));
nad_size=zeros(numel(unique_size), numel(stat_order));
rf_size=zeros(numel(unique_size), numel(stat_order));
time_size=zeros(numel(unique_size), numel(stat_order));

dup_accuracy= accuracy(dup, njtree.TruePhylo_dup);
lost_accuracy= accuracy(lost, njtree.TruePhylo_lost);
nad_accuracy= accuracy(nad, njtree.TruePhylo_nad);

for i=1:numel(unique_size)
    ind_i=(treesize==unique_size(i));
    dup_size(i, :)= (sum(dup_accuracy(ind_i, :), 1)*100.0)./nnz(ind_i);
    lost_size(i, :)= (sum(lost_accuracy(ind_i, :), 1)*100.0)./nnz(ind_i);
    nad_size(i, :)= (sum(nad_accuracy(ind_i, :), 1)*100.0)./nnz(ind_i);
    rf_size(i, :) = sum(rf(ind_i, :)==0, 1)*100.0./nnz(ind_i);
    time_size(i, :)= mean(time(ind_i, :), 1);
end

figure,
plot(unique_size, dup_size);
ylim([0, 120])
legend(stat_order);
xlabel('number of genes per tree')
title('Accuracy of inferred apparent duplication for gene trees of increasing size for simulated fungal dataset')
ylabel('apparent duplication accuracy (%)');

figure, 
plot(unique_size, lost_size);
ylim([0, 120])
legend(stat_order);
xlabel('number of genes per tree')
title('Accuracy of inferred lost for gene trees of increasing size for simulated fungal dataset')
ylabel('non apparent lost accuracy (%)');

figure, 
plot(unique_size, nad_size);
ylim([0, 120])
legend(stat_order);
xlabel('number of genes per tree')
ylabel('non apparent duplication accuracy (%)');

title('Accuracy of inferred non-apparent duplication for gene trees of increasing size for simulated fungal dataset')


figure, 
plot(unique_size, rf_size);
ylim([0, 120])
legend(stat_order);
xlabel('number of genes per tree');
ylabel('Topologies correct (%)');
title('Accuracy of the inferred topology for gene trees of increasing size for simulated fungal dataset')


figure, 
plot(unique_size, time_size(:,2:end));
ylim([0, max(max(time_size))*2])
legend(stat_order(2:end));
xlabel('number of genes per tree');
ylabel('time (s)');
title('Runtime for gene trees of increasing size for simulated fungal dataset')

end


function acc=accuracy(matrix, trueval)
    acc=bsxfun(@eq,matrix,trueval);
end



