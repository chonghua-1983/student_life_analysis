function [alpha, beta, gamma, phi] = parameter_selection(X,S,num_clu)
%num_clu = length(unique(label));
[W,H] = nndsvd(X,num_clu,0); 

knn = 12; 
flag = 1;
[L1,Dv,~] = computeHGraph_knn(X,knn,num_clu,flag);
D = (L1+abs(L1))/2; 
A = (abs(L1)-L1)/2;

n = size(S,1);
error_1 = norm(X-W*H,'fro')^2;
error_2 = norm(S-W*W','fro')^2;
error_3 = trace(W'*(D-A)*W);
error_4 = norm(W'*W- eye(num_clu), 'fro')^2;
error_5 = norm(S*ones(n,1)-ones(n,1), 'fro')^2;

alpha = error_1/error_2; 
beta = error_1/error_3;
gamma = error_1/error_4; 
phi = error_1/error_5;

% scaling parameter
% alpha = alpha;
% beta = beta*2;
% gamma = gamma;

end