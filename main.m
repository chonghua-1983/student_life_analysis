% a implementation for AIE research
% step 1: read term-doc matrix, word-occurance matrix
% term = termdoc(2:end, 1);
% term1 = wordoccurance(:,1);
% isequal(term, term1);  % 1

%% step2 
load('psychology.mat')
data_occur = termdoc1(2:end,2:end);
term_occur = termoccurance(:,2:end);
data_occur = double(data_occur); term_occur = double(term_occur);

S = ochiia(term_occur);
X = data_occur;
% Diff = dist2(X,X); 
% S = affinityMatrix(Diff);
% cluster by SNMF
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

