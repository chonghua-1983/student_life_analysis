function displayClusters(W,group,flag)
%%This function displays the clustering resultsing on the graph.
%%W is the similarity graph
%%group is the clusteirng results

%%%substract the self-similarity
if nargin<3
    flag=1;
end
W = (W +W')/2;
W = W - diag(diag(W));
%%%Normalize the similarities
n = length(W);
W(W<median(W(:)))=0;
U = unique(group);
allindex = [];
position = 0;
for i = 1 : length(U)
    index = find(group==U(i));
    allindex = [allindex;index];
    groupnew(position(i)+1:position(i)+length(index)) = index;
    position(i+1) = position(i)+length(index);
end
figure;
h = imagesc(W(allindex,allindex));hold on;
colorbar('Xtick',[1,60], 'FontSize', 12);hold on; %[1,60], 'XtickLabels',[0,1], 
if flag
    for i = 2:length(position)-1
        plot(repmat(position(i),length(W),1),1:length(W),'-w'); hold on; %,'Linewidth',3
        plot(1:length(W),repmat(position(i),length(W),1),'-w'); hold on; % ,'Linewidth',3
    end
end

h = xlabel('terms'); set(h,'FontSize',12);

h = ylabel('terms'); set(h,'FontSize',12);
h = title('Visualization on S');

set(gca, 'FontSize',12)
