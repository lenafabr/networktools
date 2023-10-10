function [NT, newregbound, boundnodeind, boundedgeind] = breakNetworkRegion(NT0,regbound)
% break network edges at the boundary of a preselected region
%regbound = nx2 coordinates defining polygonal region
% NT = new network
% newregbound = new regional boundary, including broken points
% boundnodeind = indices of new network nodes corresponding to new boundary
% points
% boundedgeind = indices of new network edges corresponding to new boundary
% edges

regsize = size(regbound,1);

if (norm(regbound(1,:)-regbound(end,:))>eps*10)
    % wrap around region boundary
    regboundw = [regbound; regbound(1,:)];
else
    regboundw = regbound;
end


NT = copy(NT0);
NT.setCumEdgeLen(1:NT.nedge,true);

% old2newind gives index to map from old regbound points to new ones
newregbound = regbound;
old2newind = 1:regsize;

boundedgeind = zeros(1,regsize);
boundnodeind = zeros(1,regsize);
for ec = 1:NT.nedge
    n1 = NT.edgenodes(ec,1); n2 = NT.edgenodes(ec,2);
    endpos = NT.nodepos([n1,n2],:); 
    
    inp = inpolygon(endpos(:,1),endpos(:,2),regbound(:,1),regbound(:,2));
    
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
    [xi,yi,ii] = polyxpoly(path(:,1),path(:,2),regbound(:,1),regbound(:,2));
    
    breakfrac = (NT.cumedgelen{ec}(ii(1)) + norm([xi yi] - path(ii(1),:)))/NT.cumedgelen{ec}(end);    
    NT.breakEdge(ec,breakfrac,1);
    NT.setCumEdgeLen(1:NT.nedge,true);
    
    % get regional boundaries with new intersections intercalated in    
    iin = old2newind(ii(2));
%     if (ii(2) == size(regbound,1))
%         newregbound = [newregbound(1:iin,:); xi yi];            
%     else
%         newregbound = [newregbound(1:iin,:); xi yi; newregbound(iin+1:end,:)];                
%     end
%           
    % which network node and edge corresponds to the boundary intersection    
    if inp(1)
        ecn = NT.nedge;       
    else
        ecn = ec;       
    end
    
    if (ii(2) == size(regbound,1))
        newregbound = [newregbound(1:iin,:); xi yi];      
        boundedgeind = [boundedgeind(1:iin) ecn];
        boundnodeind = [boundnodeind(1:iin) NT.nnode];
    else
        boundedgeind = [boundedgeind(1:iin) ecn boundedgeind(iin+1:end)];
        boundnodeind = [boundnodeind(1:iin) NT.nnode boundnodeind(iin+1:end)];
        newregbound = [newregbound(1:iin,:); xi yi; newregbound(iin+1:end,:)];              
    end
    
     old2newind(ii(2)+1:regsize) = old2newind(ii(2)+1:regsize)+1;    
 
end

end