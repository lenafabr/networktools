function [newnodepos,newedgenodes,mapold2new,mapnew2oldedge] = truncateNetworkNodes(keepind,nodepos,edgenodes)
% truncate network to only keep the node indices given by keepind
% remove any edges not between nodes in the new network.
% keep only largest connected component
% mapnew2old(i) = index of new node i in the old node list
% mapold2new(i) = index of old node i in new node list

% mapnew2oldedge(i) = for a given new edge index, gives the corresponding
% old edge index
if (isempty(keepind))
    newnodepos = [];
    newedgenodes = [];
    mapold2new = NaN*ones(size(nodepos,1),1);
    mapnew2oldedge = [];
    return
end

%%
mapold2new = NaN*ones(size(nodepos,1),1);
for ic = 1:length(keepind)
    mapold2new(keepind(ic)) = ic;
end

newnodepos = nodepos(keepind,:);
newedgenodes = [];
ct = 0;
for ec = 1:size(edgenodes,1)
    
    b1 = find(keepind==edgenodes(ec,1),1);
    b2 = find(keepind==edgenodes(ec,2),1);
    
    
    if (~isempty(b1) & ~isempty(b2))
        ct = ct+1;
        newedgenodes = [newedgenodes; b1 b2];
        mapnew2oldedge(ct) = ec;
    end
end

end