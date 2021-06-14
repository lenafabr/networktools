function NT = makeDeciHexNetwork(N, celldiam, options)
    % make a decimated honeycomb network with given parameteres
    % N, number lattice boxes
    % celldiam, cell diameter (in um)
    % options.rmnodefrac = remove certain fraction of nodes while stying
    % fully connected
    % options.rmedgefrac = remove certain fraction of edges while staying
    % fully connected
    
    opt = struct();
    opt.ntrimterminal = 0; % how many times to trim terminal nodes
    opt.rmnodefrac = 0.1;
    opt.rmedgefrac = 0.1;
    opt.dodisplay = 1 % plot at the end
    % ntrimterminal, trim this many times to remove terminal nodes
    % rmnodefrac, node removed fraction
    % rmedgefrac, edge removed fraction
    
    %copy options if given
    if (exist('options','var'))
        opt = copyStruct(options,opt);
    end
    
    % make Hex Network
    NT = makeHexNetwork(N);

    %center and scale up
    NT.nodepos = NT.nodepos-[0.5,0.5]; 
    NT.nodepos = NT.nodepos*celldiam;
    
    %remove nodes outside of a circular area
    cent = [0,0];
    diffs = NT.nodepos - cent;
    dists =sqrt(sum(diffs.^2,2));    
    keepind = find(dists<=celldiam/2);
    NT.keepNodes(keepind);
    
    %remove random nodes
    if (opt.rmnodefrac>0)
        noderm = round(opt.rmnodefrac*NT.nnode);
        [NT,whichnoderemoved]=decimateNodes(NT,noderm);
    end
    
    %decimate network edges
    if (opt.rmedgefrac>0)
        edgerm = round(opt.rmedgefrac*NT.nedge);
        [NT,whichremoved]=decimateEdges(NT,edgerm);
    end
       
    %remove terminal nodes several times
    for tc = 1:opt.ntrimterminal
        keepind = find(NT.degrees>1);
        NT.keepNodes(keepind);
    end
    
    %interpolate edgepaths
    NT.interpolateEdgePaths(2);
    % resent lengths and cumulative lengths
    NT.setCumEdgeLen(1:NT.nedge,true)
    
    %
    if (opt.dodisplay)
        NT.plotNetwork()
    end
end

