function NT = filterActiveNetwork(nodepos,edgenodes,nodeact,edgeact,edgepath)
% from a full list of nodes and edges (some inactive)
% filter out only the active ones and make a clean network object

NT = NetworkObj()

% copy over active node positions
NT.nodepos = nodepos(nodeact,:);

% copy over active edges
newedgenodes = edgenodes(edgeact,:);

% Clean up which nodes the edges are pointing to:

% this maps the list of all nodes to their index in the list of active nodes
mapall2act = zeros(1,size(nodepos,1));
mapall2act(nodeact) = 1:nnz(nodeact);

NT.edgenodes(:,1) = mapall2act(newedgenodes(:,1))
NT.edgenodes(:,2) = mapall2act(newedgenodes(:,2))

% copy over active edge paths
NT.edgepath = edgepath(edgeact);

% set up all the other arrays in the network
NT.setupNetwork();

% deal with cumulative and total edge lengths
NT.setCumEdgeLen(1:NT.nedge,true);

end