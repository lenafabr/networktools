function [NTnew,whichremoved]=decimateEdges(NT,nrm)
% starting from a network object
% attempt to remove a certain number of edges, while maintaining a single
% connected component
% returns new network object and list of original edge indices that were
% removed

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
nedge = size(G.Edges,1);

% what original index corresponds to each index in the new graph
mapnew2orig = 1:nedge;
nremoved = 0;
for rc = 1:nrm
   % [rc nrm]
    edgeorder = randperm(nedge); % order in which to try removing    
    for tc = 1:nedge % try every possible edge
        % decide which edge to remove
        erm = edgeorder(tc);
        %[rc erm]
        Gtry = G.rmedge(erm);
        [~,binsizes] = conncomp(Gtry);
        ncomp = length(binsizes); % number of connected components
        if (ncomp==1); break; end
    end
    if (ncomp>1)
        fprintf('Failed to disconnect edge %d.\n',rc)
        break
    else
        G = Gtry;
        nedge = nedge-1;
        whichremoved(rc) = mapnew2orig(erm);
        mapnew2orig(erm) = [];
        nremoved = nremoved + 1;
    end    
end

whichremoved = whichremoved(1:nremoved);
%% convert from graph back to network
NTnew = copy(NT);
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