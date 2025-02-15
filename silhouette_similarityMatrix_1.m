function [silh, s]=silhouette_similarityMatrix_1(affinity_matrix,idx,cluster)

num = length(affinity_matrix);
similar = affinity_matrix-diag(diag(affinity_matrix));
%dissimilar = 1-affinity_matrix;
%dissimilar = sqrt(1-affinity_matrix);
%dissimilar = 1-(affinity_matrix-diag(diag(affinity_matrix))+eye(num)); 
g = zeros(num,cluster); %  g(i,:)= the dissimilarity of i-th point with points in other clusters
k = zeros(1,cluster);% number of points in each cluster
for i=1:cluster
    id = find(idx==i);
    g(:,i) = sum(similar(:,id),2)/length(id);% accumulation of each row
%    k(i)=length(id)
end
a = zeros(num,1);
b1 = zeros(num,cluster-1);
for i=1:num
%    t = g(i,:)./k;
    t = g(i,:);
    a(i) = t(idx(i));
    t(idx(i))=[];
    b1(i,:)=t;
end
b = max(b1,[],2);
s = (a'-b')./abs(max([a';b']));
tmp = isnan(s); 
s(tmp) = [];
silh = mean(s);
    
    
    

    


