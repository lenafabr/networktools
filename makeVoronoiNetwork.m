function NT = makeVoronoiNetwork(latsize,rectangle,randfrac,starthex)
%% make a network from a Voronoi tesselation
% using points originally placed on a square lattice then randomly
% perturbed
% latsize = initial lattice size (approximate)
% rectangle = [lower x, lower y, width, height]
% randfrac = random perturbation (normally distributed), st dev as a
% fraction of lattice length


randshift = randfrac*latsize; % amount of random shift, as fraction of lattice size

if (~exist('starthex','var'))
    starthex =false;
end

rectpoly = [rectangle(1) rectangle(2); rectangle(1)+rectangle(3) rectangle(2);...
    rectangle(1)+rectangle(3) rectangle(2)+rectangle(4); rectangle(1) rectangle(2)+rectangle(4); ...
    rectangle(1) rectangle(2)];



if (starthex)
    %% start with hexagonal lattice
    
    R = sqrt((rectangle(3)/2)^2 + (rectangle(4)/2)^2);
    hexwidth = 2*latsize*cos(30*pi/180);
    nhex = round(2*R / hexwidth)
    %%
    NThex = makeHexNetwork(nhex);
    NThex.nodepos = NThex.nodepos.*2*R;
    %%
    % centering
    NThex.nodepos(:,1) = NThex.nodepos(:,1)-R+rectangle(1)+rectangle(3)/2;
    NThex.nodepos(:,2) = NThex.nodepos(:,2)-R+rectangle(2)+rectangle(4)/2;
    NThex.plotNetwork();
    
    %% only keep nodes in the rectangle
    keepind = find(inpolygon(NThex.nodepos(:,1),NThex.nodepos(:,2),rectpoly(:,1),rectpoly(:,2)));
    NThex.keepNodes(keepind);
%     NThex.plotNetwork()
    
    %% make random perturbations
    NThex.nodepos = NThex.nodepos + randn(NThex.nnode,2)*randshift;
    
    pts = NThex.nodepos;
elseif (startrandom) % start with random set of points within some square boundary
    mean(rectangle(3:4))/latsize
else
    % start with square lattice
    nptx = ceil(rectangle(3)/latsize);
    valsX = linspace(rectangle(1),rectangle(1)+rectangle(3),nptx);
    npty = ceil(rectangle(4)/latsize);
    valsY = linspace(rectangle(2),rectangle(2)+rectangle(4),npty);
    [X,Y] = meshgrid(valsX,valsY);
    
    pts = [X(:) Y(:)] + randn(nptx*npty,2)*randshift;
end
%% do voronoi tesselation
[v,vcells] = voronoin(pts);
rectpoly = [rectangle(1) rectangle(2); rectangle(1)+rectangle(3) rectangle(2);...
    rectangle(1)+rectangle(3) rectangle(2)+rectangle(4); rectangle(1) rectangle(2)+rectangle(4); ...
    rectangle(1) rectangle(2)];
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

end
