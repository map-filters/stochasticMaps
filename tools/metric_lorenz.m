function distMat = metric_lorenz(dimState)
%Pairwise distance matrix
distMat = zeros(dimState, dimState, 'int64');
for ii = 2:dimState
    for jj=1:ii-1
        distMat(ii,jj) =  distanceCircle(ii, jj, dimState);
    end
end
distMat = distMat + distMat';
end

function distval = distanceCircle(node1, node2, dimState)
dist1 = abs(node2-node1);
if node1<node2
   dist2 = node1+dimState - node2; 
elseif node1>node2
   dist2 = node2+dimState - node1;  
else
    dist2=0;
end
distval = min(dist1, dist2);
end