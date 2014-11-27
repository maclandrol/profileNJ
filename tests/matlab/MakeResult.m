datatype='fungi';
param='all';
fungi11=load('fungi11');
fungi22=load('fungi22');
fungi41=load('fungi41');
fungi44=load('fungi44');
fungi11rac=load('fungi11rac');
fungi22rac=load('fungi22rac');
fungi41rac=load('fungi41rac');
fungi44rac=load('fungi44rac');
%load([datatype, param, '.mat'])
%load([datatype, param,'rac', '.mat'])
%%concat at the same time

nj100=vertcat(fungi11.nj100, fungi22.nj100, fungi41.nj100, fungi44.nj100);
nj101=vertcat(fungi11.nj101, fungi22.nj101, fungi41.nj101, fungi44.nj101);
nj90=vertcat(fungi11.nj90, fungi22.nj90, fungi41.nj90, fungi44.nj90);
nj95=vertcat(fungi11.nj95, fungi22.nj95, fungi41.nj95, fungi44.nj95);
nj80=vertcat(fungi11.nj80, fungi22.nj80, fungi41.nj80, fungi44.nj80);
nj75=vertcat(fungi11.nj75, fungi22.nj75, fungi41.nj75, fungi44.nj75);
nj50=vertcat(fungi11.nj50, fungi22.nj50, fungi41.nj50, fungi44.nj50);
njlkl=vertcat(fungi11.njlkl, fungi22.njlkl, fungi41.njlkl, fungi44.njlkl);
njtree=vertcat(fungi11.njtree, fungi22.njtree, fungi41.njtree, fungi44.njtree);

%nracine data
njractree=vertcat(fungi11rac.njractree, fungi22rac.njractree, fungi41rac.njractree, fungi44rac.njractree);
njraclkl=vertcat(fungi11rac.njraclkl, fungi22rac.njraclkl, fungi41rac.njraclkl, fungi44rac.njraclkl);
njrac95=vertcat(fungi11rac.njrac95, fungi22rac.njrac95, fungi41rac.njrac95, fungi44rac.njrac95);


%Star script

true_dlc=[njtree.TruePhylo_dup+njtree.TruePhylo_nad,njtree.TruePhylo_lost, njtree.TruePhylo_dlc];

%all matrix data
all_ad_matrix=[njtree.Raxml_dup, njtree.Treefix_dup, nj50.PolySolver_dup, nj75.PolySolver_dup,  nj90.PolySolver_dup, nj95.PolySolver_dup, nj100.PolySolver_dup, nj101.PolySolver_dup, njrac95.PolySolver_dup];
all_nad_matrix=[njtree.Raxml_nad, njtree.Treefix_nad, nj50.PolySolver_nad, nj75.PolySolver_nad,  nj90.PolySolver_nad, nj95.PolySolver_nad, nj100.PolySolver_nad, nj101.PolySolver_nad, njrac95.PolySolver_nad];
all_lost_matrix=[njtree.Raxml_lost, njtree.Treefix_lost, nj50.PolySolver_lost, nj75.PolySolver_lost,  nj90.PolySolver_lost, nj95.PolySolver_lost, nj100.PolySolver_lost, nj101.PolySolver_lost, njrac95.PolySolver_lost];
all_recon_matrix=[njtree.Raxml_dlc, njtree.Treefix_dlc, nj50.PolySolver_dlc, nj75.PolySolver_dlc,  nj90.PolySolver_dlc, nj95.PolySolver_dlc, nj100.PolySolver_dlc, nj101.PolySolver_dlc, njrac95.PolySolver_dlc];
all_order={'RAxML', 'TreeFix','profileNJ-50', 'profileNJ-75', 'profileNJ-90', 'profileNJ-95', 'profileNJ-100', 'profileNJ-*', 'profileNJ rooted'};
all_rf_matrix=[njtree.Raxml_rf, njtree.Treefix_rf, nj50.PolySolver_rf, nj75.PolySolver_rf,  nj90.PolySolver_rf, nj95.PolySolver_rf, nj100.PolySolver_rf, nj101.PolySolver_rf, njrac95.PolySolver_rf];
all_maxrf_matrix=[njtree.Raxml_maxrf, njtree.Treefix_maxrf, nj50.PolySolver_maxrf, nj75.PolySolver_maxrf,  nj90.PolySolver_maxrf, nj95.PolySolver_maxrf, nj100.PolySolver_maxrf, nj101.PolySolver_maxrf, njrac95.PolySolver_maxrf];

%cellfun( @(x) sprintf(x), all_order, 'UniformOutput', false)

% data to show
ad_matrix=[njtree.Raxml_dup,njtree.Treefix_dup, nj95.PolySolver_dup];
nad_matrix=[njtree.Raxml_nad,njtree.Treefix_nad, nj95.PolySolver_nad];
lost_matrix=[njtree.Raxml_lost,njtree.Treefix_lost, nj95.PolySolver_lost];
recon_matrix= [njtree.Raxml_dlc,njtree.Treefix_dlc, nj95.PolySolver_dlc];
rf_matrix = [njtree.Raxml_rf, njtree.Treefix_rf, nj95.PolySolver_rf];
maxrf_matrix = [njtree.Raxml_maxrf, njtree.Treefix_maxrf, nj95.PolySolver_maxrf];
lkl_matrix=[njlkl.Raxml_logL, njlkl.TreeFix_logL, njlkl.ProfileNJ95_logL, njlkl.TruePhylo_logL];
au_matrix=[njlkl.Raxml_AU, njlkl.TreeFix_AU, njlkl.ProfileNJ95_AU, njlkl.TruePhylo_AU];
mltime=nj95.ML_time;
racmltime=njrac95.ML_time;
mltime(nj95.polysolver_nsol<2)=0;
racmltime(njrac95.polysolver_nsol<2)=0;

time= [njtree.Treefix_time, njrac95.PolySolver_time+racmltime, nj95.PolySolver_time, nj95.PolySolver_time+mltime];
order={'RAxML', 'TreeFix','profileNJ', 'Expected Tree'};

time_label={'TreeFix', 'profileNJ root', 'profileNJ Unroot', 'profileNJ + logL comp.'};

%% Phase 1 : recon cost
reconcost(true_dlc, all_order, all_ad_matrix+all_nad_matrix, all_lost_matrix, all_recon_matrix);
reconcost(true_dlc, order(1:end-1), ad_matrix+nad_matrix, lost_matrix, recon_matrix);

%% Phase 2 : Rf computation
RFfigure(rf_matrix, order(1:end-1));
    % RF for profileNJ alpha comparision
RFfigure(all_rf_matrix(:, 3:end) , all_order(3:end));

%% Phase 3: rf by cost computation
%rfbycost(all_rf_matrix, all_maxrf_matrix, true_dlc(:,end), all_order, 5)
rfbycost(rf_matrix, maxrf_matrix, true_dlc(:,end), order(1:end-1), 2)

%% Phase 4: likelihood
likelihood(lkl_matrix, au_matrix, rf_matrix, maxrf_matrix, order)

%% Phase 5 : Size stat
sizeStat(njtree.leaves, njtree, ad_matrix,lost_matrix, nad_matrix,  time, rf_matrix, order(1:end-1), time_label, datatype, 10);


%% Time xx nsol
timesol([nj95.PolySolver_time, nj95.PolySolver_time+mltime], nj95.polysolver_nsol, {'profileNJ', 'profileNJ + logL comp.'},2);

savefig(fullfile('images', strcat(datatype, param)), strcat(datatype, param))

close all;


% tabulate is best rf
tabulate(nj95.is_bestrf)
