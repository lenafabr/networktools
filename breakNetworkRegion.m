function [NT, newregbound, newboundnodeind, mapold2new] = breakNetworkRegion(NT0,regbound,options)
% break network edges at the boundary of a preselected region
%regbound = nx2 coordinates defining polygonal region
% NT = new network
% newregbound = new regional boundary, including broken points
% boundnodeind = indices of new network nodes corresponding to new boundary
% points
% boundedgeind = indices of new network edges corresponding to new boundary
% edges


% default options
opt = struct();
opt.keepin = true; % keep nodes within network?
opt.keepout = true; % keep nodes outside network?
% do not make any new edges smaller than this length
% will adjust boundary to avoid this where necessary
opt.minedgelen = 0;
% only allow degree 1 nodes to be connected to broken region
opt.condeg1only = false;

if (exist('options','var'))
    opt = copyStruct(options,opt);
end

regsize = size(regbound,1);

if (norm(regbound(1,:)-regbound(end,:))>eps*10)
    % wrap around region boundary
    regboundw = [regbound; regbound(1,:)];
else
    regboundw = regbound;
end


NT = copy(NT0);
NT.setCumEdgeLen(1:NT.nedge,true);

nodesinregion = inpolygon(NT.nodepos(:,1),NT.nodepos(:,2),regbound(:,1),regbound(:,2));

% old2newind gives index to map from old regbound points to new ones
newregbound = regbound;
old2newind = 1:regsize;

boundedgeind = zeros(1,regsize);
boundnodeind = zeros(1,regsize);
boundparam0 = arclenparam(regboundw);
boundparam = boundparam0;

% keep these nodes even if they are not in the desired region
stillkeep = false(NT.nnode,1);

for ec = 1:NT.nedge
    n1 = NT.edgenodes(ec,1); n2 = NT.edgenodes(ec,2);
    endpos = NT.nodepos([n1,n2],:); 
    
    %inp = inpolygon(endpos(:,1),endpos(:,2),regbound(:,1),regbound(:,2));
    inp = nodesinregion([n1 n2]);
    
    if (inp(1) & ~inp(2))
        innode = n1; outnode = n2;
    elseif (inp(2) & ~inp(1))
        innode = n2; outnode = n1;
    else
        continue;
    end
       
    %plot(NT.nodepos(innode,1),NT.nodepos(innode,2),'r*')
    
    % get path intersection
    path = NT.edgepath{ec};    
    [xi,yi,ii] = polyxpoly(path(:,1),path(:,2),regboundw(:,1),regboundw(:,2));
    
    % parameter of intersection along the boundary
    newboundparam = boundparam0(ii(2)) + norm([xi,yi] - regboundw(ii(2),:));

    elen = NT.cumedgelen{ec}(end);
    breakfrac = (NT.cumedgelen{ec}(ii(1)) + norm([xi yi] - path(ii(1),:)))/elen;

    if (elen < opt.minedgelen)
        % edge is already too short. Do not break further        
        dobreak = false;

        if (breakfrac<0.5)
            newnode = n1;
        else
            newnode = n2;
        end
        if (opt.condeg1only & NT.degrees(newnode)>1)
            error('PROBLEM: cannot connect sheet at node of degree >1. Manually adjust sheet')
        end

        boundpt = NT.nodepos(newnode,:);
        stillkeep(newnode) = true;
    else
        dobreak = true;
        
        if (~opt.keepout)
            if (inp(1) & breakfrac*elen<opt.minedgelen)
                breakfrac = opt.minedgelen/elen;
            elseif (inp(2) & (1-breakfrac)*elen<opt.minedgelen)
                breakfrac = 1 - opt.minedgelen/elen;
            end
        end
        if (~opt.keepin)
            if (inp(2) & breakfrac*elen<opt.minedgelen)
                breakfrac = opt.minedgelen/elen;
            elseif (inp(1) & (1-breakfrac)*elen<opt.minedgelen)
                breakfrac = 1 - opt.minedgelen/elen;
            end
        end

        % break the edge
        NT.breakEdge(ec,breakfrac,1);
        NT.setCumEdgeLen(1:NT.nedge,true);

        newnode = NT.nnode;
        % position of new point
        boundpt = NT.nodepos(end,:);
    end

    % get regional boundaries with new intersections intercalated in    
    iin = old2newind(ii(2));
%     if (ii(2) == size(regbound,1))
%         newregbound = [newregbound(1:iin,:); xi yi];            
%     else
%         newregbound = [newregbound(1:iin,:); xi yi; newregbound(iin+1:end,:)];                
%     end
%           
    % which network node and edge corresponds to the boundary intersection    
    if (inp(1) & dobreak)
        ecn = NT.nedge;       
    else
        ecn = ec;       
    end
    
    if (ii(2) == size(regbound,1))
        newregbound = [newregbound(1:iin,:); boundpt];      
        boundedgeind = [boundedgeind(1:iin) ecn];
        boundnodeind = [boundnodeind(1:iin) newnode];
        boundparam = [boundparam(1:iin); newboundparam];
    else
        boundedgeind = [boundedgeind(1:iin) ecn boundedgeind(iin+1:end)];
        boundnodeind = [boundnodeind(1:iin) newnode boundnodeind(iin+1:end)];
        newregbound = [newregbound(1:iin,:); boundpt; newregbound(iin+1:end,:)];    
        boundparam = [boundparam(1:iin); newboundparam; boundparam(iin+1:end)];
    end
    
     old2newind(ii(2)+1:regsize) = old2newind(ii(2)+1:regsize)+1;     
end

% sort the region boundary points in appropriate order (arclength around
% original region)
[boundparam,bi] = sort(boundparam(1:end-1));
boundnodeind = boundnodeind(bi);
boundedgeind = boundedgeind(bi);
newregbound = newregbound(bi,:);

% remove nodes outside region? inside region?
keepind = 1:NT.nnode;

if (~opt.keepout | ~opt.keepin)
    if (~opt.keepout)
        keepind(find(~nodesinregion & ~stillkeep)) = [];
    elseif (~opt.keepin)
        keepind(find(nodesinregion & ~stillkeep)) = [];
    end

    mapold2newedge = zeros(1,NT.nedge);
    [mapold2new,mapnew2oldedge] = NT.keepNodes(keepind);
    
    ind = find(boundnodeind>0);
    newboundnodeind = zeros(size(boundnodeind));
    newboundnodeind(ind) = mapold2new(boundnodeind(ind));

    ind = find(mapnew2oldedge>0);
    mapold2newedge(mapnew2oldedge(ind)) = ind;

    ind = find(boundedgeind>0);
    newboundedgeind = zeros(size(boundedgeind));
    newboundedgeind(ind)=  mapold2newedge(boundedgeind(ind));
else
    newboundnodeind = boundnodeind;
    newboundedgeind = boundedgeind;
end

end