function NTnew = joinNetworks(NT1,NT2,newedges)
% join together two networks
% newedges = [i,j]: additional edges between node i in NT1 and node j in NT2

NTnew = NetworkObj();

if (~exist('newedges','var'))
    newedges = []; % no new edges by default
end

if (NT1.dim~=NT2.dim)
    error('cannot join two networks with different dimensions')
end

nn1 = NT1.nnode;

NTnew.nodepos = [NT1.nodepos; NT2.nodepos];

% update indices in original edgenodes
edges2 = NT2.edgenodes+nn1;
% add new edges
if (~isempty(newedges))
    newedges(:,2) = newedges(:,2)+nn1;
end

NTnew.edgenodes = [NT1.edgenodes; edges2; newedges];

% set up new network
NTnew.setupNetwork()


end