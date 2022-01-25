load('../examples/exampleERnetwork.mat')

%% Create a networkgraph object 'from scratch'
NTg = NetworkGraphObj();
NTg.setupNetwork(NT.edgenodes,NT.nodepos,NT.edgepath');
NTg0 = copy(NTg);
%NTg.keepNodes([50,35,54,55])

%%
NTg.redistributeEdgePaths(4);

%%
spNT = adjustableNetwork(NTg);

%% convert back to old network to check
NTold = NTg.convert2OldNetwork();

%% play with gui
app = networkGraphEdit('nt',NTg,'img',img)

%%
appchildren = get(app.NetworkGraphEditGUIUIFigure,'Children');
for c = 1:length(appchildren)
    set(appchildren(c),'Enable','off')
end


%% try adding edges
paths{1} = [1 2; 3 4];
paths{2} = [5 6; 7 8];
NTg.addEdges([203;249],[227;206])

%% try merging at nodes
NTg = copy(NTg0);
NTg.mergeEdgesAtNode(243)
NTg0.plotNetwork(struct('nodecolor',[1 0 0]));
hold all
NTg.plotNetwork()
hold off
set(gca,'YDir','reverse')
%% play with splineroi
points = [1 2; 3 5; 2 6; 1 7];
spobj = splineroi();
spobj.addNode(points)

%% makea small network
NTg.keepNodes([132,152,172,180,179,157,151])

%% play with spline network
spNT = splineNetwork(NTg)

%%
roipts = findobj(app.plotfig,'Type','images.roi.point');
goodind = cellfun(@(x) ~isempty(x), {roipts.Position})
roipts = roipts(goodind);

for rc = 1:length(roipts)
    pos = roipts(rc).Position;
    hold all
    plot(pos(1),pos(2),'m*','MarkerSize',30)    
end
hold off

%% try selecting edges
app.graphplotH.PickableParts = 'none';
for ec = 1:length(app.edgeplotH)
    app.edgeplotH(ec).PickableParts= 'all';
end


%% Set up new type of network object

load('../examples/exampleERnetwork.mat')

%% Create a networkgraph object 'from scratch'
NTg = NetworkGraphObj();
NTg.setupNetwork(NT.edgenodes,NT.nodepos,NT.edgepath');
NTg0 = copy(NTg);
%% plot graph
plotopt = struct('showEdgeDir',0,'edgelabels',1:NTg.graph.numedges);
[graphplotH,edgeplotH] = NTg.plotNetwork(plotopt);

%% 
NTg = copy(NTg0);
nodepos = [10,10;11, 11];
connect{1} = [1,2,3];
connect{2} = [4,5,6];
%%
NTg.addNodes(nodepos,connect)

%%
NTg.breakEdge(25,30,0.5)

% ----------- OLD STUFF -----------

%% try adding edge paths to graphs
G = NT.makeGraph();
G.Edges.edgepath = NT.edgepath';
G.Nodes.XData = NT.nodepos(:,1);
G.Nodes.YData = NT.nodepos(:,2);

%%

for ec = 1:NT.nedge
    edgepath = NT.edgepath{ec};
    edgeh(ec) = plot(edgepath(:,1),edgepath(:,2),'-');
    hold all
end
graphh = G.plot('XData',G.Nodes.XData,'YData',G.Nodes.YData,'EdgeAlpha',0)
labeledge(graphh,1:NT.nedge,1:NT.nedge)
hold off

%% Try making new networkgraph object
NTg = NetworkGraphObj(G)