%% example script to make hexagonal (honecomb) networks 
% in a circular domain, with central circle cut out to represent nucleus
% optionally, can truncate some nodes or edges

celldiam = 40; % cell diameter (in um)
nucdiam = 15; % nucleus diameter (in um)

N =20; % number lattice boxes along one dimension (density of honeycomb)

% make honeycomb network
NT0 = makeHexNetwork(N);

% visualize it
NT0.plotNetwork()
title('Original honeycomb network')

% ------
%% expand to cell perimeter and cut out nucleus

NT = copy(NT0);

% center and scale up
NT.nodepos = NT.nodepos-[0.5,0.5];
NT.nodepos = NT.nodepos*celldiam;
dists = sqrt(sum(NT.nodepos.^2,2));

keepind = find(dists<=celldiam/2 & dists>=nucdiam/2);

% keep only nodes in the annular region, and the edges that connect them
NT.keepNodes(keepind);

% visualize remaining network
NT.plotNetwork()
title('Annular honeycomb network')

%% if desired, get rid of some nodes and edges to make a less well-connected network
nodermfrac = 0.05; % remove this fraction of the nodes
edgermfrac = 0.05; % remove this fraction of the edges

% remove network nodes, while attempting to maintain single connected
% component (might remove fewer, if it can't find a way to maintain
% connectivity)
noderm = round(nodermfrac*NT.nnode);
[NT,whichnoderemoved]=decimateNodes(NT,noderm);

% remove network edges, while attempting to maintain single connected
% component (might remove fewer, if it can't find a way to maintain
% connectivity)
edgerm = edgermfrac*NT.nedge;
[NT,whichremoved]=decimateEdges(NT,edgerm);

% visualize remaining network
NT.plotNetwork()
title('Decimated honeycomb network')

%% save network to file
NT.outputNetwork('example_periphhex.net')