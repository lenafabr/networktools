% Build a simple triskelion network object as an example

%% define points in space
nodepos = [0 0; 0 1; -1 -1; 2 -1];
%% edge nodes list pairs of nodes that are connected to each other
edgenodes = [1 2; 1 3; 1 4];

%% build a network object
NT = NetworkObj()
NT.nodepos = nodepos;
NT.edgenodes = edgenodes;
NT.setupNetwork()

NT.plotNetwork(struct('labels',1))

%% mark a node as permeable
NT.nodelabels{1} = 'P1';

%% output to a network object file
NT.outputNetwork('triskelion_example.net')