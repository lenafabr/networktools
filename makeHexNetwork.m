function NT = makeHexNetwork(N)
% build a hexagonal lattice
% N = number of honeycomb units along each dimension

%%
r = 1; % center to corner of hex

w = r*cos(pi/6); % half-width of hex

% first row of points
y = 0
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
for hc = 1:N-1
    newpts = bsxfun(@plus,firstpts,[0,hc*(3*r)]);
    allpts = [allpts; newpts];    
end
% last row
newpts = bsxfun(@plus,pt1,[0,N*(2*r+w)]);
allpts = [allpts; newpts];

% keep only points within square
%R = N*2*w;
%allpts = allpts/R;
ind = find(allpts(:,1)<N*2*w & allpts(:,2)<N*2*w*0.99);
allpts = allpts(ind,:);

%plot(allpts(:,1),allpts(:,2),'.')
%axis equal

%% rescale to width 1
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
%NT.plotNetwork()