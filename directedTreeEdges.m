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

            % adjust edge widths if necessary
            % the 2nd column of edgewidth contains length along edge where
            % the width measurement is made
            widths = NT.edgewidth{ec};
            if ~isempty(widths)
                if (size(widths,2)>1)
                    NT.edgewidth{ec}(:,2) = NT.edgelens(ec) - widths(:,2);
                end
            end
        else
            error('problem with edge node connectivity')
        end   
        
        % fix nodeedges so that the parent edge is always first
        nd = NT.edgenodes(ec,2);
        tmpedgelist = NT.nodeedges(nd,:);
        i = find(tmpedgelist==ec);
        tmpedgelist(i) = [];                
        NT.nodeedges(nd,:) = [ec tmpedgelist];
        tmpnodelist = NT.nodenodes(nd,:);
        tmpnodelist(i) = [];
        NT.nodenodes(nd,:) = [nc tmpnodelist];       
        
        % call recursively starting from the next node
        [isset,wasreversed] = directedTreeEdges(NT,nextnode,isset,wasreversed);
    end
end

    
end