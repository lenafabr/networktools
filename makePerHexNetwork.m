function [NT,newedge] = makePerHexNetwork(N)
% build a periodic hexagonal lattice
% N = number of honeycomb units along each dimension
% scaled so that horizontal and vertical size is roughly the same

r = 1; % center to corner of hex

w = r*cos(pi/6); % half-width of hex

% first row of points
y = 0;
pt1 = [(0:N)'*2*w y*ones(N+1,1)];

% 2nd row of points
pt2 = bsxfun(@plus,pt1,[w,r*sin(pi/6)]);

% 3rd row of points
pt3 = bsxfun(@plus,pt2,[0 r]);

% 4th row of pts
pt4 = bsxfun(@plus,pt1,[0,2*r]);

firstpts = [pt1; pt2; pt3; pt4];

% repeat rows several times
allpts = firstpts;
for hc = 1:(round(sqrt(3)*N/3))
    newpts = bsxfun(@plus,firstpts,[0,hc*(3*r)]);
    allpts = [allpts; newpts];    
end
% last row
% newpts = bsxfun(@plus,pt1,[0,N*(2*r+w)]);
% allpts = [allpts; newpts];

% keep only points within square
%R = N*2*w;
%allpts = allpts/R;
% ind = find(allpts(:,1)<N*2*w & allpts(:,2)<N*2*w*0.99);
% allpts = allpts(ind,:);

% plot(allpts(:,1),allpts(:,2),'o')
% axis equal

% rescale to width 1
maxx = max(allpts(:,1));
allpts = allpts/maxx;

NT = NetworkObj();
NT.nodepos = allpts;
NT.nnode = size(NT.nodepos,1);

% connect up nearest neighbors
ct=0;
for pc1 = 1:NT.nnode
    for pc2 = pc1+1:NT.nnode
        if (norm(NT.nodepos(pc1,:)-NT.nodepos(pc2,:))< r/maxx+1e-8)
            ct = ct+1;
            NT.edgenodes(ct,:) = [pc1,pc2];
        end
    end
end

NT.nedge = size(NT.edgenodes,1);
NT.setupNetwork()

% NT.plotNetwork()

% edgevals tracks which edges are periodicity connections
NT.edgevals = zeros(NT.nedge,1);
%% connect boundary points together:
nedge = NT.nedge;
% first row to final row (bottom to top)
maxy = max(NT.nodepos(:,2));

nodesRow1 = find(NT.nodepos(:,2)==0);
nodesRowN = find(NT.nodepos(:,2)>=maxy-1e-8);
% connect them
for pc = 1:length(nodesRow1)
    NT.edgenodes(end+1,:) = [nodesRow1(pc) nodesRowN(pc)];
    NT.edgelens(end+1) = NT.edgelens(1); % edgelengths will be the same as all other edges
    NT.edgepath{end+1} = [NT.nodepos(nodesRow1(pc),:) ; NT.nodepos(nodesRowN(pc),:)];
    ec = size(NT.edgenodes,1);
    NT.edgevals(ec) = 2;
end
NT.nedge=size(NT.edgenodes,1);

%% find first and final columns
maxx=max(NT.nodepos(:,1));
nodesCol1 = find(NT.nodepos(:,1)==0);
nodesColN = find(NT.nodepos(:,1)>=maxx-1e-8);

for pc=1:length(nodesCol1)
    NT.edgenodes(end+1,:) = [nodesCol1(pc) nodesColN(pc)];
    NT.edgelens(end+1) = NT.edgelens(1); % edgelengths will be the same as all other edges
    NT.edgepath{end+1} = [NT.nodepos(nodesCol1(pc),:) ; NT.nodepos(nodesColN(pc),:)];
    ec = size(NT.edgenodes,1);
    NT.edgevals(ec) = 1;
end
NT.nedge=size(NT.edgenodes,1);
NT.edgevals(end+1:NT.nedge) = 0;
NT.setupNetwork % don't reset edgelengths, but update degrees, anything else that may have been missed


% how many new edges appended at the end. These are real edges, but their
% spatial dependence is wonky, have to be careful with them
newedge=NT.nedge-nedge;

end