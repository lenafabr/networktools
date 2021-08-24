function [NT,skelimage,opt] = getNetworkFromBWImage(bwimg,options)


% skeletonize a B&W image and extract a network structure
% input: 
% bwimg = black and white image (nxm array of integer 0s and 1s)
% options = structure containing any optional parameters
% output: 
% NT = network structure
% skelimage = skeletonized image;

%% default options
opt = struct();
% how much displaying to do during the calculation
opt.dodisplay = 0;

% if set to true, keep only the largest connected component
opt.keepconncomp = true;

% maximum allowed jump in an edge, in pixels
opt.maxedgestep = 5; 

% trim any terminal edges shorter than this (in pixels)
opt.mintermlen = 3;

% how close a node has to be to break a contour
opt.closenode = 2;

% copy input parameters
if (exist('options','var'))
    opt = copyStruct(options,opt);
end
%% get skeletonized image
skelimage=bwmorph(bwimg,'skel',inf);
   
%% get branch nodes and end nodes
branchimg = bwmorph(skelimage,'branchpoints');
endimg = bwmorph(skelimage,'endpoints');
nodeimg = branchimg + endimg;
nodeimg(nodeimg>1) = 1;

clear nodepoints
[nodepoints(:,2),nodepoints(:,1)] = ind2sub(size(skelimage),find(nodeimg));

%%
if (opt.dodisplay>1) % show nodes on skeleton image
    figure
    imshowpair(skelimage,nodeimg)
    hold all
    plot(nodepoints(:,1),nodepoints(:,2),'co')
    hold off
    drawnow
end
    
%% edges
edgeimg = skelimage-nodeimg;
edgeimg(edgeimg<0) = 0;

comp = bwconncomp(edgeimg);

edgesize = cellfun(@(x) length(x),comp.PixelIdxList);


%% group directly connected nodes together
compnodes = bwconncomp(nodeimg,4);

for cc = 1:compnodes.NumObjects
    pix = compnodes.PixelIdxList{cc};
    if (length(pix)==1)
        [nodegroup(cc,2),nodegroup(cc,1)] = ind2sub(size(skelimage),pix);
    else
        [y,x] = ind2sub(size(skelimage),pix);
        nodegroup(cc,1) = mean(x);
        nodegroup(cc,2) = mean(y);
    end
end

%% get rid of redundant nodes
duplicatenodes = false(size(nodegroup,1),1);
for ncc = 1:size(nodegroup,1)-1
    diffs = nodegroup(ncc+1:end,:)-nodegroup(ncc,:);
    dist = sum(diffs.^2,2);
    badind = find(dist==0);
    duplicatenodes(ncc+badind) = true;
end
nodegroup(duplicatenodes,:) =[];
%%

%
if (opt.dodisplay>1) % show grouped nodes on image
    figure
    imshowpair(skelimage,nodeimg)
    hold all
    plot(nodegroup(:,1),nodegroup(:,2),'co')   
    hold off
    drawnow
end
    
%% starting from a node, trace out all connected edges
% directions to step in
step = [1 1; 1 0; 1 -1; 0 1; 0 -1; -1 1; -1 0; -1 -1];
nd = size(step,1);

checkmat = zeros(size(skelimage));
edgect = 0;
alledges= {};
edgenodes = [];

remainimg= edgeimg;
for nc = 1:size(nodegroup,1)
    %[nc size(nodegroup,1)]
    %%
    subs = [round(nodegroup(nc,2))+step(:,2),round(nodegroup(nc,1))+step(:,1)];
    
    %
    good = subs(:,1)>0 & subs(:,1)<=size(skelimage,1) & subs(:,2)>0 & subs(:,2) <=size(skelimage,2);
    subs = subs(good,:);
    inds = sub2ind(size(skelimage),subs(:,1),subs(:,2));
    nodeind =  sub2ind(size(skelimage),round(nodegroup(nc,2)),round(nodegroup(nc,1)));
    
    % directions we can go in
    tryind = find(skelimage(inds));
    dirs = inds(tryind);

    dirlist = {'SE','E','NE','S','N','SW','W','NW'};
    dirlistuse = dirlist(tryind);
    %% trace each edge
    tmpimg = remainimg;
    tmpimg(nodeind) = 1;
    
    for dc= 1:length(dirs)
        %%
        useimg = tmpimg;
        useimg(dirs) =0;
        useimg(dirs(dc)) = 1;
        
        
        contour = bwtraceboundary(useimg,fliplr(round(nodegroup(nc,:))),dirlistuse{dc});
        % avoid doubling back
        endind = round((size(contour,1)+1)/2);
        contour = contour(1:endind,:);
        
        stopind = size(contour,1);
        for cc = 1:size(contour,1)
            %%
            dists = sqrt(sum((contour(cc,:) - fliplr(nodegroup)).^2,2));
            dists(nc) = inf;
            closenodes = find(dists<opt.closenode);
            if (~isempty(closenodes))
                % break contour
                nearestnode = closenodes(1);
                stopind = cc;
                contour(cc,:) = fliplr(nodegroup(nearestnode,:));
                break
            end
            
            if (cc==size(contour,1))
                [~,nearestnode] = min(dists);
            end
        end
        contour = contour(1:stopind,:);
        contour(1,:) = fliplr(nodegroup(nc,:));
        
        edgect= edgect+1;
        alledges{edgect} = [fliplr(contour); nodegroup(nearestnode,:)];
        edgenodes(edgect,:) = [nc nearestnode];
        
        % knockout edges from image
        rmind = sub2ind(size(skelimage),contour(2:end-1,1),contour(2:end-1,2));
        remainimg(rmind)= 0;
        
    end
    
end
%toc

%% show nodes and edges
if (opt.dodisplay>1)
    figure
    imshow(skelimage)
    hold all
    plot(nodegroup(:,1),nodegroup(:,2),'co')
    cmap = lines(edgect);
    for ec = 1:edgect
        edge = alledges{ec};
        plot(edge(:,1),edge(:,2),'.-','MarkerSize',10,'Color',cmap(ec,:))
        nodes = edgenodes(ec,:);
        plot(nodegroup(nodes,1),nodegroup(nodes,2),'r.','MarkerSize',20)
        text(edge(1,1),edge(1,2),sprintf('%d',ec),'Color',cmap(ec,:))
    end
    hold off
    drawnow
end
%% sort node order for each edge
for ec = 1:edgect
    n1 = edgenodes(ec,1); n2 = edgenodes(ec,2);
    if (n2<n1)
        edgenodes(ec,:) = [n2 n1];
        alledges{ec} = flipud(alledges{ec});
    end
end

%% clean up jumpy redundant edges
nnode = size(nodegroup,1);
ct = 0;
saveedges = {};
saveedgenodes = [];

for n1 = 1:nnode
    % [n1 nnode]
    useedge = (edgenodes(:,1)==n1);
    for n2 = n1+1:nnode
        useedge2  = useedge & (edgenodes(:,2)==n2);
        
        useind = find(useedge2);
        if (~isempty(useind))
            maxstepsize = zeros(length(useind),1);
            for cc = 1:length(useind)
                edge = alledges{useind(cc)};
                diffs = diff(edge);
                stepsize = sqrt(sum(diffs.^2,2));
                maxstepsize(cc) = max(stepsize);
            end
            
            % keep connection with lowest stepsize
            % and only if its below a cutoff
            [minval,savecc] = min(maxstepsize);
            if (minval<opt.maxedgestep)
                ct = ct+1;
                saveedges{ct} = alledges{useind(savecc)};
                saveedgenodes(ct,:)= [n1, n2];
            end
        end
    end
end


%% remove backtracks along edges
cleanedges = saveedges;
for ec = 1:length(saveedges)
    while 1
        edge = cleanedges{ec};
        newedge = edge;
        didcontract = false;
        for cc = 1:size(edge,1)
            diffs = edge(cc+1:end,:)-edge(cc,:);
            dists = sum(diffs.^2,2);
            
            repeatind = find(dists==0,1);
            
            if (~isempty(repeatind))
                newedge(cc+1:cc+repeatind,:) = [];
                didcontract = true;
                break
            end
        end
        cleanedges{ec} = newedge;
        if (~didcontract)
            break
        end
    end
end
cleanedgenodes = saveedgenodes;


%% show nodes and edges
if (opt.dodisplay>1) 
    figure
    imshow(skelimage)
    hold all
    plot(nodegroup(:,1),nodegroup(:,2),'co')
    cmap = lines(length(cleanedges));
    for ec = 1:length(cleanedges)
        edge = cleanedges{ec};
        plot(edge(:,1),edge(:,2),'.-','MarkerSize',10,'Color',cmap(ec,:))
        nodes = edgenodes(ec,:);
        plot(nodegroup(nodes,1),nodegroup(nodes,2),'r.','MarkerSize',20)
        text(edge(1,1),edge(1,2),sprintf('%d',ec),'Color',cmap(ec,:))
    end
    hold off
    drawnow
end
    
%% Smooth edges using smoothing splines
edgestep = 1; % step size for defining curvy edges

if (opt.dodisplay>0) % show final network
    cmap = lines(length(cleanedges));
    figure 
    imshow(skelimage)
    hold all
    plot(nodegroup(:,1),nodegroup(:,2),'go')
end
clear smoothedges
for ec = 1:length(cleanedges)
    edge = cleanedges{ec};
    param = arclenparam(edge');
    nstep = ceil(param(end)/edgestep);
    newparam = linspace(param(1),param(end),nstep);
    
    %npt = ceil(length/edgestep)
    pp = csaps(param,edge',0.5);
    
    smedge = ppval(pp,newparam)';
    % plot(edge(:,1),edge(:,2),'g.-')
    if(opt.dodisplay>0)
        plot(smedge(:,1),smedge(:,2),'.-','Color',cmap(ec,:),'LineWidth',2)
    end
    smoothedges{ec} = smedge;
end
if (opt.dodisplay>0)
    hold off
    drawnow
end

%% make a network object
NT = NetworkObj();
NT.dim = 2;
NT.nodepos = nodegroup;
NT.edgenodes = cleanedgenodes;
NT.edgepath = smoothedges;

NT.setupNetwork();
NT.removeDoubleEdges();
NT.setupNetwork();
NT.setEdgeLens();
NT.edgepath = smoothedges;
NT.setCumEdgeLen();

%% trim short terminal nodes

dokeep = false(1,NT.nnode);
for nc = 1:NT.nnode
    if NT.degrees(nc)==1
        ec = NT.nodeedges(nc,1);
        if (NT.edgelens(ec)>opt.mintermlen)
            dokeep(nc) = 1;
        end
    else
        dokeep(nc) = 1;
    end
end

keepind = find(dokeep);
NT.keepNodes(keepind);
NT.setupNetwork()

if (opt.keepconncomp)
    NT.keepLargestConnComp()   
end

% adjust edge paths to actually hit nodes
for ec = 1:NT.nedge
    NT.edgepath{ec}(1,:) = NT.nodepos(NT.edgenodes(ec,1),:);
    NT.edgepath{ec}(end,:) = NT.nodepos(NT.edgenodes(ec,2),:);
end

NT.setCumEdgeLen(1:NT.nedge,1);

%% empty edgewidths arrays for use later
NT.edgewidth = cell(NT.nedge,1);

