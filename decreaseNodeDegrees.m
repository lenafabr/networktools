function decreaseNodeDegrees(NT,shiftfrac,maxct)
% shift edges on a network so that there are no nodes of degree above 3
% for each high degree node, randomly choose one edge
% and shift its endpoint off that node to another connected edge
% shifted by fraction shiftfrac of that edge's length

% maxct = maximum number of edges to fix

highdegnodes = find(NT.degrees>3);
ct = 0;

if (~exist('maxct','var'))
    maxct = inf;
end

while ~isempty(highdegnodes)
    ct = ct+1;
       
    
    if (ct>maxct)
        disp('stop shifting edges')
        break
    end
    
    nhigh = length(highdegnodes);
    nc = highdegnodes(randi(nhigh));
    
    if (mod(ct,100)==0); disp([ct nhigh]); end
    
    %[nc nhigh]
    
    deg = NT.degrees(nc);
    % pick which edge to shift and which to connect it to
    tmp = randsample(NT.nodeedges(nc,1:deg),2);
    ecshift = tmp(1); ectarget = tmp(2);
    
    % create a new node along target edge, breaking up that edge
    if (NT.edgenodes(ectarget,1)==nc) % outgoing edge
        NT.breakEdge(ectarget,shiftfrac,0);
    elseif (NT.edgenodes(ectarget,2)==nc) % incoming edge
        NT.breakEdge(ectarget,1-shiftfrac,0);
    else
        error('connectivity problem')
    end
    
    % shift edge
    if (NT.edgenodes(ecshift,1)==nc) % outgoing edge
        NT.edgenodes(ecshift,1) = NT.nnode+1;
    elseif (NT.edgenodes(ecshift,2)==nc) % incoming edge
        NT.edgenodes(ecshift,2) = NT.nnode+1;
    else
        error('connectivity problem')
    end
    
    NT.setupNetwork();
    alterededges = [ecshift,ectarget,NT.nedge];
    NT.interpolateEdgePaths(2,alterededges);
    NT.setCumEdgeLen(alterededges,1);
    
    highdegnodes = find(NT.degrees>3);
end

end