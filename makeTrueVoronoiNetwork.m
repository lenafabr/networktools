function NT = makeTrueVoronoiNetwork(Npts)

%%
x = rand([1 Npts]);
y = rand([1 Npts]);
pts=[x' y'];

%% do voronoi tesselation
[v,vcells] = voronoin(pts);
rectpoly = [0 0 ; 1 1];
inrect = inpolygon(v(:,1),v(:,2),rectpoly(:,1),rectpoly(:,2));
% 
% plot(pts(:,1),pts(:,2),'b.')
% hold all
% plot(v(insquare,1),v(insquare,2),'r.')
% hold off

% get edges from cells
A = zeros(size(v,1)); % adjacency matrix
for cc = 1:length(vcells)
    inds = vcells{cc};
    for ec = 1:length(inds)-1
        A(inds(ec),inds(ec+1)) = 1;
        A(inds(ec+1),inds(ec)) = 1;
    end
    A(inds(end),inds(1)) = 1;
    A(inds(1), inds(end)) = 1;
end

ct = 0;
for ic = 1:size(v,1)
    for jc = ic+1:size(v,1)        
        if (A(ic,jc) & inrect(ic) & inrect(jc))
            ct=ct+1;
            edgenodes(ct,:) = [ic jc];
        end
    end
end

%% make network
NT = NetworkObj();
NT.nodepos = v;
NT.edgenodes = edgenodes;
NT.setupNetwork()
keepind = find(inrect);
NT.keepNodes(keepind);
NT.setupNetwork()
NT.setEdgeLens();

% plot(pts(:,1),pts(:,2),'r.')
% hold all
% NT.plotNetwork()
% hold off
%%
end