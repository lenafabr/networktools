% code for making a stretched out hex network roughly representing ATL KO
% structures

%% load in the ATL network regions for comparison
dirname='../networks/ATLKO/';
load([dirname,'ATLKO_circreg.mat'])

dirname='../networks/WT/';
load([dirname,'WT_circreg.mat'])

% --------------
%% get WT edge lengths
%% load in WT networks
dirname='../networks/WT/';
load([dirname,'WT_circreg.mat'])
dirname='../networks/ATLKO/';
load([dirname,'ATLKO_circreg.mat'])


%%
edgelenWT = [];
for fc = 1:length(WTnetworks)
    NT = copy(WTnetworks(fc));
    NT.mergeAllEdgePaths(); % merge away degree 2 nodes
    NT.setCumEdgeLen(1:NT.nedge,true);
    % pick out only edges between junctions
    %ind = find(NT.degrees(NT.edgenodes(:,1))> 2 & NT.degrees(NT.edgenodes(:,2))>2);
    edgelenWT = [edgelenWT; NT.edgelens];
end
avglenWT = mean(edgelenWT)

%%
edgelenATL = [];
for fc = 1:length(ATLnetworks)
    NT = copy(ATLnetworks(fc));
    NT.mergeAllEdgePaths(); % merge away degree 2 nodes
    NT.setCumEdgeLen(1:NT.nedge,true);
    %ind = find(NT.degrees(NT.edgenodes(:,1))> 2 & NT.degrees(NT.edgenodes(:,2))>2);
    edgelenATL = [edgelenATL; NT.edgelens];
end
avglenATL = mean(edgelenATL)



%% load in honeycomb network
NT0 = NetworkObj('../networks/circlehexR25N25.net',struct('dim',2));

% rescale to avg edge len in WT network
NT0.scaleCoords(avglenWT/mean(NT0.edgelens));
NT = copy(NT0);
%%
NT0.plotNetwork()

%% stretch out the network to an appropriate aspect ratio and rough polygon size
NTstr = copy(NT0);
%NTstr.nodepos = NTstr.nodepos.*[1.5,3];
NTstr.nodepos = NTstr.nodepos.*[1,2.32];
NTstr.interpolateEdgePaths(2);
NTstr.setCumEdgeLen(1:NTstr.nedge,true);

%NTstr.scaleCoords(avglenATL/mean(NTstr.edgelens))
%%
NTstr.plotNetwork()

%% Cut out a circular region
th = linspace(0,2*pi,50)';
nodedist = sum(NTstr.nodepos.^2,2);
[~,centnode] = min(nodedist);
centpos = NTstr.nodepos(centnode,:)-[0,1];
cellR = 8; % cell radius
circbound = centpos + cellR*[cos(th) sin(th)];

NTstr.plotNetwork()
hold all
plot(circbound(:,1),circbound(:,2),'g.-')
hold off
%%
[NT, newregbound, boundnodeind, boundedgeind] = breakNetworkRegion(NTstr,circbound,struct('keepout',false,'minedgelen',0.3));

%% trim degree 1 nodes on short edges
dokeep = true(1,NT.nnode);
for nc = find(NT.degrees==1)'
    ec = NT.nodeedges(nc,1);
    elen = NT.edgelens(ec);

    if (elen < 1);
        dokeep(nc) = false;
    end
end
NT.keepNodes(find(dokeep));

%% make a special release node
nodedist = sum(NT.nodepos.^2,2);
[~,centnode] = min(nodedist);
n2 = NT.nodenodes(centnode,1);
breakpos = (NT.nodepos(centnode,:)+NT.nodepos(n2,:))/2;
[distmin, breakfracmin] = setNodeNearestPoint(NT,breakpos,0.3);

%%
NT.plotNetwork(struct('datatipindex',true))
hold all
plot(circbound(:,1),circbound(:,2),'g.-')

plot(NT.nodepos(end,1),NT.nodepos(end,2),'r*')
hold off

%% set node labels and output
NT.nodelabels = {};
for nc = 1:NT.nnode
    NT.nodelabels{nc} = '';
end
NT.nodelabels{end} = 'P1';

NT.outputNetwork('../networks/circlehexstretch2pt3_elen1pt4.net')

% -----------
%% Make a corresponding size network without the stretch

th = linspace(0,2*pi,50)';
nodedist = sum(NT0.nodepos.^2,2);
[~,centnode] = min(nodedist);
centpos = NT0.nodepos(centnode,:);
cellR = 8; % cell radius
circbound = centpos + cellR*[cos(th) sin(th)];

[NT, newregbound, boundnodeind, boundedgeind] = breakNetworkRegion(NT0,circbound,struct('keepout',false,'minedgelen',0.3));

%% trim degree 1 nodes on short edges
dokeep = true(1,NT.nnode);
for nc = find(NT.degrees==1)'
    ec = NT.nodeedges(nc,1);
    elen = NT.edgelens(ec);

    if (elen < 0.5);
        dokeep(nc) = false;
    end
end
NT.keepNodes(find(dokeep));

%% make a special release node
nodedist = sum(NT.nodepos.^2,2);
[~,centnode] = min(nodedist);
n2 = NT.nodenodes(centnode,1);
breakpos = (NT.nodepos(centnode,:)+NT.nodepos(n2,:))/2;
[distmin, breakfracmin] = setNodeNearestPoint(NT,breakpos,0.3);

%%
NT.plotNetwork(struct('datatipindex',true))
hold all
plot(circbound(:,1),circbound(:,2),'g.-')

plot(NT.nodepos(end,1),NT.nodepos(end,2),'r*')
hold off

%% set node labels and output
NT.nodelabels = {};
for nc = 1:NT.nnode
    NT.nodelabels{nc} = '';
end
NT.nodelabels{end} = 'P1';

NT.outputNetwork('../networks/circlehexR8deg2P_elenpt83.net')

%%
save('../networks/circlehexstretch2pt3_elen.mat')
