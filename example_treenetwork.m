%% Setting up an example tree network

NT = NetworkObj();

% set up node positions
NT.nodepos = [0 2; 0 1; -1 0; 1 0];

% set up connectivity
NT.edgenodes = [
    2 1
    3 2
    4 2
    ];

%% set up network
NT.setupNetwork()
%% get edge lens
NT.setEdgeLens(true)

%% plot the network
NT.plotNetwork(struct('labels',true));

%% list some info
% node degrees:
degrees = NT.degrees

% edge lengths
edgelens = NT.edgelens

% which edges connect to each node?
nodeedges = NT.nodeedges

