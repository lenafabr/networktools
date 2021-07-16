function [NTnew,whichremoved]=decimateNodes(NT,nrm,keepnodes)
% starting from a network object
% attempt to remove a certain number of nodes, while maintaining a single
% connected component
% returns new network object and list of original edge indices that were
% removed

if (~exist('keepnodes','var'))
    keepnodes = [];
end

%% convert network object to matlab graph object

% get adjacency matrix
A = zeros(NT.nnode,NT.nnode);
for ec = 1:NT.nedge
    nodes = NT.edgenodes(ec,:);
    A(nodes(1),nodes(2)) = 1;
    A(nodes(2),nodes(1)) = 1;
end

G = graph(A);
[bins,binsizes] = conncomp(G);
ncomp = length(binsizes); % number of connected components

if (ncomp>1)
    error('original network is already disconnected')
end

%%
nnode = size(G.Nodes,1);

% what original index corresponds to each index in the new graph
mapnew2orig = 1:nnode;
nremoved = 0;

canremove = true(1,nnode);
canremove(keepnodes)  = false;

for rc = 1:nrm
   % [rc nrm]
    nodeorder = randperm(nnode); % order in which to try removing    
    success = false;
    for tc = 1:nnode % try every possible node
        % decide which node to remove
        erm = nodeorder(tc);
        if (~canremove(mapnew2orig(erm)))
            continue
        end
        %[rc erm]
        Gtry = G.rmnode(erm);
        [~,binsizes] = conncomp(Gtry);
        ncomp = length(binsizes); % number of connected components
        if (ncomp==1)
            success = true;
            break
        end
    end
    if (~success)
        fprintf('Failed to disconnect node %d.\n',rc)
        break
    else
        G = Gtry;
        nnode = nnode-1;
        whichremoved(rc) = mapnew2orig(erm);
        mapnew2orig(erm) = [];
        nremoved = nremoved + 1;
    end    
end

whichremoved = whichremoved(1:nremoved);
%% convert from graph back to network
NTnew = NT;
NTnew.nnode = nnode;
NTnew.nodepos(whichremoved,:) = [];
NTnew.edgenodes = G.Edges.EndNodes;
NTnew.edgelens(whichremoved) = [];
if (~isempty(NTnew.edgepath))
    NTnew.edgepath(whichremoved) = [];
end
if (~isempty(NT.edgevals))
    NTnew.edgevals(whichremoved) = [];
end
NTnew.setupNetwork(false);

end