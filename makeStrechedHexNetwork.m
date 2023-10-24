% code for making a stretched out hex network roughly representing ATL KO
% structures

%% load in the ATL network regions for comparison
dirname='../ERCaSims/networks/ATLKO/';
load([dirname,'ATLKO_circreg.mat'])

dirname='../ERCaSims/networks/WT/';
load([dirname,'WT_circreg.mat'])

%%
for fc = 1:length(ATLnetworks)
    avgedgelen(fc) = mean(ATLnetworks(fc).edgelens);
end
for fc = 1:length(WTnetworks)
    avgedgelenWT(fc) = mean(WTnetworks(fc).edgelens);
end

%% load in honeycomb network
NT0 = NetworkObj('../ERCaSims/networks/circlehexR25N25.net',struct('dim',2));
NT = copy(NT0);
%%
NT0.plotNetwork()

%% stretch out the network to an appropriate aspect ratio and rough polygon size
NTstr = copy(NT0);
NTstr.nodepos = NTstr.nodepos.*[1.5,3];
NTstr.interpolateEdgePaths(2);
NTstr.setCumEdgeLen(1:NTstr.nedge,true);
%%
NTstr.plotNetwork()

%% Cut out a circular region
th = linspace(0,2*pi,50)';
nodedist = sum(NTstr.nodepos.^2,2);
[~,centnode] = min(nodedist);
centpos = NTstr.nodepos(centnode,:);
cellR = 15; % cell radius
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
for nc = 1:NT.nnode
    NT.nodelabels{nc} = '';
end
NT.nodelabels{end} = 'P1';

NT.outputNetwork('../ERCaSims/networks/circlehexstretch3.net')

% -----------
%% Make a corresponding size network without the stretch

th = linspace(0,2*pi,50)';
nodedist = sum(NT0.nodepos.^2,2);
[~,centnode] = min(nodedist);
centpos = NT0.nodepos(centnode,:);
cellR = 15; % cell radius
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
for nc = 1:NT.nnode
    NT.nodelabels{nc} = '';
end
NT.nodelabels{end} = 'P1';

NT.outputNetwork('../ERCaSims/networks/circlehexR15deg2P.net')
