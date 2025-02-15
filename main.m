% a implementation for AIE research
% step 1: read term-doc matrix, word-occurance matrix
term = termdoc(2:end, 1);
% term1 = wordoccurance(:,1);
% isequal(term, term1);  % 1

%% step2 
load('psydata/psychology.mat')
data_occur = termdoc1(2:end,2:end);
term_occur = termoccurance(:,2:end);
data_occur = double(data_occur); term_occur = double(term_occur);

S = ochiia(term_occur);
X = data_occur;
% Diff = dist2(X,X); 
% S = affinityMatrix(Diff);
% cluster by SNMF
addpath('symnmf-master')
D = diag(sum(S)); L = D-S;
[K1, K2, K12,K22] = Estimate_Number_of_Clusters_given_Laplacian(L,(2:10));

num_clu = K12; 
% [X, ~] = tfidf2(X);
knn = 12; flag = 1;
[L,Dv,Av] = computeHGraph_knn(X,knn,num_clu,flag);
HyperD = (L+abs(L))/2; hyperA = (abs(L)-L)/2;

[W, H] = nndsvd(X, num_clu, 0);
S = S./repmat(sum(S),size(S,1),1);
[alpha, beta, gamma, phi] = parameter_selection(X, S, num_clu);

[Wb,Hb,Sb,objs] = HSNMF(X, W, H, S, HyperD, hyperA, 0, beta, gamma);
A = Wtrim(Sb,20); 
[clust,~,~] = getNCluster(A, num_clu, 0,3,20);
[sil_ave_HSNMF, sil_HSNMF] = silhouette_similarityMatrix_1(A, clust, num_clu);
displayClusters(A, clust);
% ch metric
addpath('cvik-toolbox-master/proximity')
addpath('cvik-toolbox-master/cvi')
ch_HSNMF = chindex(clust,X,'DISTANCE','lap');
db_HSNMF = dbindex(clust,X,'DISTANCE','lap');

% visualization on W
addpath('nmfv1_4')
option.propertyName='threshold';
option.propertyValue=0.9;
[ACluster,YCluster,indACluster,indYCluster,Xout,Aout,Yout,tElapsed] = bicluster_hsnmf(X,Wb,Hb,option);
option.colormap=winter; % redgreencmap summer hsv hot cool spring winter autumn
option.standardize=true;
NMFBicHeatMap(Xout,Aout,Yout,indACluster,indYCluster,num_clu,option)







% tsne 
% reduction = tsne(Wb);
reduction = tsne_p(A);
term = '.'; 
colors = generateColors(max(length(unique(clust))),length(unique(term)));
gscatter(reduction(:,1),reduction(:,2),clust,colors,term,9); %_5k_10k
set(gca,'xtick',[],'ytick',[]);
title('HSNMF')
% UMAP
addpath('umapFileExchange/umap') 
addpath('umapFileExchange/util')
term = '.'; 
colors = generateColors(max(num_clu,length(unique(term))));
% clust = label_name;

D_orthnnmf = 1-A; D_orthnnmf = D_orthnnmf-diag(diag(D_orthnnmf));
[reduction_hsnmf, ~, ~, ~] = run_umap(D_orthnnmf,'metric','precomputed','min_dist',0.25,'n_neighbors',8); %,'min_dist',0.2,'n_neighbors',5
% 'min_dist',0.3,'n_neighbors',5  NMFGOT;
% NMFGOT on H
% [reduction_hsnmf, ~, ~, ~] = run_umap(Wb,'min_dist',0.78,'n_neighbors',12); %,'min_dist',0.68,'n_neighbors',12
gscatter(reduction_hsnmf(:,1),reduction_hsnmf(:,2),clust,colors,[],10);
set(gca,'xtick',[],'ytick',[]);
title('UMAP for OrthNMF')
legend('Location','westoutside','Box','off','FontSize',9.5); 
legendmarkeradjust(16)
legend('boxoff') 

%% baseline comparison
% NMF
[wbest,hbest,normbest] = nnmf(X,num_clu);
rng(2019);
idx = randperm(size(X,1));
[indic, center, ~, ~, ~] = litekmeans(wbest, num_clu, 'MaxIter', 50, 'Start', idx(1:num_clu)');
sil_NMF = silhouette(X, indic);  % 'cosine' Euclidean
sil_ave_NMF = mean(sil_NMF);
ch_NMF = chindex(indic,X,'DISTANCE','lap');
db_NMF = dbindex(indic,X,'DISTANCE','lap');

% SC
% Sim = ochiia(term_occur);
Diff = dist2(X,X); 
Sim = affinityMatrix(Diff,50);
[group, ~] = SpectralClustering(Sim, num_clu);
ch_sc = chindex(group,X,'DISTANCE','lap');
db_sc = dbindex(group,X,'DISTANCE','lap');
sil_sc = silhouette(X, group);  % 'cosine' Euclidean
sil_ave_sc = mean(sil_sc);

% DBSCAN-----not stable
% eps = 0.1; mpts = 8;
% idx = dbscan(X,eps,mpts);
% ch_dbscan = chindex(idx,X,'DISTANCE','euc');
% db_dbscan = dbindex(idx,X,'DISTANCE','euc');
% sil_dbscan = silhouette(X, idx);  % 'cosine' Euclidean
% sil_ave_dbscan = mean(sil_dbscan);

% SNMF
[idx, iter, obj, Wsy] = symnmf_cluster(X, num_clu); 
ch_sym = chindex(idx,X,'DISTANCE','lap');
db_sym = dbindex(idx,X,'DISTANCE','lap');
sil_sym = silhouette(X, idx);  % 'cosine' Euclidean
sil_ave_sym = mean(sil_sym);

%% term topics -- visualization on S
topics = cell(num_clu,1);
for i=1:num_clu
    tmp = find(clust == i);
    topics{i} = term(tmp);
end


