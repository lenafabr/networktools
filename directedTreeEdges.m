function [isset,wasreversed] = directedTreeEdges(NT,topnode,isset,wasreversed)
% change around edge directions on a network to convert to a directed tree
% starting at the given topnode
% edges point down the tree
% will also fix edgepaths that got reversed previously


if (all(isset))
    % completely done setting tree
    % recalculate cumulative edge lengths
    NT.setCumEdgeLen(1:NT.nedge,true)
end

nc = topnode; % current node

% adjacent edges
adjedges = NT.nodeedges(nc,1:NT.degrees(nc));
% adjacent edges that aren't set yet
unsetedges = adjedges(~isset(adjedges));


if isempty(unsetedges) % done with this node, everything is set
    return
else
    for ec = unsetedges
        if (NT.edgenodes(ec,1)==nc) % edge is already correct
            isset(ec) = true;
            nextnode = NT.edgenodes(ec,2);
            wasreversed(ec) = false;
        elseif (NT.edgenodes(ec,2)==nc) % edge needs flipping
            nextnode = NT.edgenodes(ec,1);
            NT.edgenodes(ec,:) = fliplr(NT.edgenodes(ec,:));
            NT.edgepath{ec} = flipud(NT.edgepath{ec});            
            isset(ec) = true;
            wasreversed(ec) = true;
        else
            error('problem with edge node connectivity')
        end   
        % call recursively starting from the next node
        [isset,wasreversed] = directedTreeEdges(NT,nextnode,isset,wasreversed);
    end
end

    
end