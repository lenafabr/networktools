% class definition for a network object
% containing data on network structure and basic utilities for
% working with networks

%classdef NetworkObj < handle
classdef NetworkObj < matlab.mixin.Copyable   
properties                               
    nnode 
    nedge 
    dim 
    degrees % degree of each node
    maxdeg % max degree in network
    nodepos % position of each node
    edgenodes % nodes connected by each edge
    nodenodes % nodes connected to each node
    nodeedges % edges connected to each node
    edgeedges % edges connected to each edge
    edgelens % edge lengths
    nodelabels % string label for each node (permeability or reservoir)
    nodevals
    edgevals
    rootnode

    % used for gui only    
    loops
    Name
    
    edgepath % coordinates of path along each edge
    cumedgelen % cumulative lengths along each edge
    edgewidth;
    
    use_edgeedge;
end

methods
    function NT = NetworkObj(fname,options)
        % create a network object
        % optionally, load data from file
        
        NT.nnode = 0; % number of nodes
        NT.nedge = 0; % number of edges
        NT.dim = 2; % spatial dimension        
        NT.maxdeg = 10; % maximum allowed degree              
        NT.loops = [];
        NT.Name = '';
        NT.rootnode = NaN;

        % keep track of edge-edge connectivity?
        NT.use_edgeedge=true;
        
        if (exist('fname','var'))
            if (exist('options','var'))
                NT.loadNetwork(fname,options);
            else
                NT.loadNetwork(fname);
            end
        end
    end
    
    function loadNetwork(NT,fname,options)
        opt = struct();
        opt.rmduplicate = 1; % remove duplicate edges
        opt.dim = 0; % dimension is predefined?
        
        if (nargin > 2)
            opt = copyStruct(options, opt);
        end
        
        % load a network from file
        [nodepos,edgenodes,edgevals,nodelabels] = loadnetworkstruct(fname);
        
        
        if (opt.dim>0)
            NT.nodepos = nodepos(:,1:opt.dim);
            NT.nodevals = nodepos(:,opt.dim+1:end);
            NT.dim = opt.dim;
        else
            NT.nodepos = nodepos;
            NT.dim = size(nodepos,2);
        end
        
        NT.edgenodes = edgenodes(:,1:2);
        NT.edgevals = edgevals;
        NT.nodelabels = nodelabels;
        
        NT.setupNetwork();
        
        % get rid of doubled edges
        if opt.rmduplicate
            NT.removeDoubleEdges();
            NT.setupNetwork() % reset arrays
        end    
    end
    
    function setupNetwork(NT,resetedgepath)
        
        if (~exist('resetedgepath','var'))
            % recalculate edge paths?
            resetedgepath = false;
        end
       
        
        % still calculate edge lens if none defined at all
        doedgelens = false; 
        
        % set up network structure based on nodepos and edgenodes arrays
        % only
        NT.dim = size(NT.nodepos,2);
        NT.nnode = size(NT.nodepos,1);
        NT.nedge = size(NT.edgenodes,1);
        
        % get degree and nodes attached to each node
        NT.degrees = zeros(NT.nnode,1);
        NT.nodenodes = zeros(NT.nnode,NT.maxdeg);
        NT.nodeedges = zeros(NT.nnode,NT.maxdeg);
        if (isempty(NT.edgelens))
            doedgelens = true;
            NT.edgelens = zeros(NT.nedge,1);
        end
                
        
        for ec = 1:NT.nedge
            n1 = NT.edgenodes(ec,1); n2 = NT.edgenodes(ec,2);
            NT.degrees(n1) = NT.degrees(n1)+1;
            NT.nodenodes(n1,NT.degrees(n1)) = n2;
            NT.nodeedges(n1,NT.degrees(n1)) = ec;
            
            NT.degrees(n2) = NT.degrees(n2)+1;
            NT.nodenodes(n2,NT.degrees(n2)) = n1;
            NT.nodeedges(n2,NT.degrees(n2)) = ec;
            
            % get edge length
            %NT.edgelens(ec) = norm(NT.nodepos(n1,:) - NT.nodepos(n2,:));
        end
        
        NT.maxdeg = max(NT.degrees);
        NT.nodenodes = NT.nodenodes(:,1:NT.maxdeg);
        NT.nodeedges = NT.nodeedges(:,1:NT.maxdeg);
        
        
        if (NT.use_edgeedge)
            % set edgeedge connections
            NT.edgeedges = zeros(NT.nedge,2,max(NT.degrees)*2);
            for ec = 1:NT.nedge
                n1 = NT.edgenodes(ec,1); n2 = NT.edgenodes(ec,2);
                % edgeedge array gives, for a given edge:
                % which node is a neighbor edge attached to (1 or 2)
                % what is the index of that neighbor edge
                ct = 0;
                for cc = 1:NT.degrees(n1)
                    ec2 = NT.nodeedges(n1,cc); % adjacent edge
                    if (ec~=ec2)
                        ct = ct+1;
                        NT.edgeedges(ec,:,ct) = [1,ec2];
                    end
                end
                for cc = 1:NT.degrees(n2)
                    ec2 = NT.nodeedges(n2,cc); % adjacent edge
                    if (ec~=ec2)
                        ct = ct+1;
                        NT.edgeedges(ec,:,ct) = [2,ec2];
                    end
                end
            end
        end
        
        if (resetedgepath || doedgelens)            
            % reset edge paths and lengths
            NT.edgepath = cell(NT.nedge,1);
            NT.edgelens = zeros(NT.nedge,1);
            for ec = 1:NT.nedge      
                n1 = NT.edgenodes(ec,1); n2 = NT.edgenodes(ec,2);
                % get default edge length; MUST be at least distance
                % between points                
                len = norm(NT.nodepos(n1,:) - NT.nodepos(n2,:));
                NT.edgelens(ec) = max(len,NT.edgelens(ec));                
            end
            NT.interpolateEdgePaths(2);
            NT.setCumEdgeLen();
        end
            
        if (isempty(NT.nodelabels))
            NT.nodelabels = cell(NT.nnode,1);
        end
               
    end
    
    function removeDoubleEdges(NT)
        %% get rid of all duplicate edges, so you only have one edge
        % between each pair of nodes
        badind = [];
        newedgenodes = [];
        edgenodes= NT.edgenodes;
        for ec1 = 1:NT.nedge
            duplicatefound = 0;
            for ec2 = ec1+1:NT.nedge
                if (edgenodes(ec1,1)==edgenodes(ec2,2) & edgenodes(ec1,2)==edgenodes(ec2,1))
                    disp([ec1 ec2 edgenodes(ec1,:) edgenodes(ec2,:)])
                    duplicatefound = 1;
                elseif (edgenodes(ec1,1)==edgenodes(ec2,1) & edgenodes(ec1,2)==edgenodes(ec2,2))
                    disp([ec1 ec2 edgenodes(ec1,:) edgenodes(ec2,:)])
                    duplicatefound=1;
                end
            end
            if ~duplicatefound
                newedgenodes = [newedgenodes; edgenodes(ec1,:)];
            end
        end
        NT.edgenodes = newedgenodes;
    end
    
    function connectNodeNearest(NT,pos)
        % connect a new node at position pos
        % to the nearest edge in the network
        
        % find nearest segment
        mindist = inf;
        for ec = 1:NT.nedge
            n1 = NT.edgenodes(ec,1); n2 = NT.edgenodes(ec,2);
            [dist(ec),frac(ec),ptnear(ec,:)] = ptlinesegdist(pos,NT.nodepos(n1,:),NT.nodepos(n2,:));
            
            if (dist(ec)<mindist)
                mindist = dist(ec);
                minec = ec;
            end
        end
        
        % create new node along segment
        n1 = NT.edgenodes(minec,1); n2 = NT.edgenodes(minec,2);
        pt1 = NT.nodepos(n1,:); pt2 = NT.nodepos(n2,:);
        newpos = pt1 + frac(minec)*(pt2-pt1);
        NT.nodepos(end+1,:) = newpos;
        NT.nnode = NT.nnode+1;
        
        % reconnect edges
        NT.edgenodes(minec,:) = [n1,NT.nnode];
        NT.edgenodes(end+1,:) = [NT.nnode,n2];
        
        % new node for target
        NT.nodepos(end+1,:) = pos;
        NT.nnode = NT.nnode+1;
        NT.edgenodes(end+1,:) = [NT.nnode-1,NT.nnode];
        
        NT.setupNetwork();
        
    end
    
    function [mapold2new,mapnew2oldedge] = keepNodes(NT,keepind)        
        % truncate network to only keep the node indices given by keepind
        % remove any edges not between nodes in the new network.
        % keep only largest connected component
        % mapnew2old(i) = index of new node i in the old node list
        % mapold2new(i) = index of old node i in new node list
        
        % for nodes: mapnew2old would be = keepind
        
        [newnodepos,newedgenodes,mapold2new,mapnew2oldedge] = truncateNetworkNodes(keepind,NT.nodepos,NT.edgenodes);
        
        NT.nodepos = newnodepos;
        NT.edgenodes = newedgenodes;
        
        NT.setupNetwork()                
                
        if (~isempty(NT.nodevals)); NT.nodevals = NT.nodevals(keepind); end
        if (~isempty(NT.edgevals)); NT.edgevals = NT.edgevals(mapnew2oldedge); end
        if (~isempty(NT.edgewidth)); NT.edgewidth = NT.edgewidth(mapnew2oldedge); end
        if (~isempty(NT.edgepath)); NT.edgepath = NT.edgepath(mapnew2oldedge); end
        if (~isempty(NT.edgelens)); NT.edgelens = NT.edgelens(mapnew2oldedge); end
        if (~isempty(NT.cumedgelen)); NT.cumedgelen = NT.cumedgelen(mapnew2oldedge); end
    end
    
    function [mapnew2oldedge] = keepEdges(NT,keepind)        
        % keep only edges of a certain index        
        if (size(keepind,1)>size(keepind,2))
            keepind = keepind';
        end

        % mapping from new to old edge index
        mapnew2oldedge = zeros(1,nnz(keepind));
        ct = 0;
        for ec = keepind
            ct = ct+1;
            mapnew2oldedge(ct) = ec;
        end
        
        NT.edgenodes = NT.edgenodes(keepind,:);        
        NT.setupNetwork()      
        
        if (~isempty(NT.edgewidth)); NT.edgewidth = NT.edgewidth(mapnew2oldedge); end               
        if (~isempty(NT.edgevals)); NT.edgevals = NT.edgevals(mapnew2oldedge); end
        if (~isempty(NT.edgepath)); NT.edgepath = NT.edgepath(mapnew2oldedge); end
        if (~isempty(NT.edgelens)); NT.edgelens = NT.edgelens(mapnew2oldedge); end
        if (~isempty(NT.cumedgelen)); NT.cumedgelen = NT.cumedgelen(mapnew2oldedge); end
    end
    
    function addNodes(NT,nodepos,connect)
        % add nodes to the network
        % connecting each to predefined old node indices
        nnew = size(nodepos,1);
        
        nnprev = NT.nnode;
        neprev = NT.nedge;
        
        nind = NT.nnode;
        eind = NT.nedge;            
        NT.nodepos = [NT.nodepos; nodepos];
                
        for pc = 1:nnew
            nind = nind+1;
            newedge = length(connect{pc});            
            for ec = 1:newedge % add in new edges to this node
                eind = eind+1;
                con = connect{pc}(ec); % which old node does it connect to?
                NT.edgenodes(eind,:) = [con, nind];
            end
        end
        
        % set up other network info arrays, based on nodepos and edgenodes
        NT.setupNetwork();        
        
        if (~isempty(NT.nodevals)); NT.nodevals(nnprev+1:NT.nnode) = 0; end
        if (~isempty(NT.edgevals)); NT.edgevals(neprev+1:NT.nedge) = 0; end        
        if (~isempty(NT.nodelabels))
            for nc = nnprev:NT.nnode
                NT.nodelabels{nc} = '';                
            end
        end
        if (~isempty(NT.edgewidth))
            for ec= neprev:NT.nedge
                NT.edgewidth{ec} = [];
            end
        end        
        
    end
    
    function outputNetwork(NT,outfile,options)
        % output network structure to file
        % edgepaths: lists of cycles, in terms of edges to include as LOOP
        % structures
        % reservoirnodes: which nodes belong to a certain reservoir
        
        opt = struct();
        opt.nodelabels = {}; % are there any labels for the nodes?
        opt.WRITEEVS = false; % should the edge values be written?
         % output extra lines for edge path coordinates?
        opt.WRITEPATHS = false;

        if (nargin > 2)
            opt = copyStruct(options, opt);
        end
        
        of = fopen(outfile,'w')
        fprintf(of,'%s\n','# file defining network structure')
        fprintf(of,'%s\n\n','# made with matlab NetworkObj output')
        
        fprintf(of,'%s \n','# list of node indices and xy positions and values')        
        nodefmtstring = ['NODE %d ' repmat(['%20.10f '],1,NT.dim+1) '\n'];
        nodelblfmtstring = ['NODE %d ' repmat(['%20.10f '],1,NT.dim+1) '%s \n'];
        %
        for pc = 1:size(NT.nodepos,1)
            
            if (isempty(NT.nodevals))
                val = 0.0;              
            else
                val = NT.nodevals(pc);                
            end
            
            if (~isempty(opt.nodelabels))
                fprintf(of,nodelblfmtstring,pc, NT.nodepos(pc,:), val,opt.nodelabels{pc});
            elseif (~isempty(NT.nodelabels))
                fprintf(of,nodelblfmtstring,pc, NT.nodepos(pc,:), val,NT.nodelabels{pc});
            else
                fprintf(of,nodefmtstring,pc, NT.nodepos(pc,:), val);
            end
        end
        % edge information
        nev = size(NT.edgevals,2);
        if(nev > 0 && opt.WRITEEVS == true)
            edgefmtstring = "EDGE %d %d %d %20.10f";
            for ev = 1:nev
                edgefmtstring = edgefmtstring+" %20.10f";
            end
            edgefmtstring = edgefmtstring+"\n";
            for ec = 1:size(NT.edgenodes)
                fprintf(of,edgefmtstring,ec, NT.edgenodes(ec,:), NT.edgelens(ec), NT.edgevals(ec,:));
            end
        else
            edgefmtstring = "EDGE %d %d %d %20.10f\n";
            for ec = 1:size(NT.edgenodes,1)
                fprintf(of,edgefmtstring,ec, NT.edgenodes(ec,:), NT.edgelens(ec));
            end
        end

        if (opt.WRITEPATHS)
            % output lines giving cartesian coordinates for each edge
            % format:
            % EDGEPATH EC NPT x1, x2, ... x_npt, y1, y2, ... y_npt, [z_1,
            % z_2, ... z_npt]
            % where NPT is the number of points along the edge, z only
            % included if network dimensionality is 3

            for ec = 1:NT.nedge
                path = NT.edgepath{ec};
                npt = size(path,1);
                fmtstr = ['EDGEPATH %d %d ' repmat(['%20.10f '],1,npt*NT.dim) '\n'];
                fprintf(of,fmtstr,ec, npt, path(:,1)',path(:,2)');
            end
        end
        
        % independent cycles as a list of edges
%         if (exist('edgepaths','var'))
%             for pc = 1:length(edgepaths)
%                 fprintf(of,['LOOP ' sprintf('%d ',edgepaths{pc}) '\n'])
%             end   
%         elseif (~isempty(NT.loops))
%             for pc = 1:length(NT.loops)
%                 fprintf(of,['LOOP ' sprintf('%d ',NT.loops{pc}) '\n'])
%             end   
%         end
        fclose(of)
    end
    
    function keepLargestConnComp(NT)
        % truncate the network, keeping only the largest connected
        % component
        
        % get graph structure
        [NTgraph,A] = NT.makeGraph();
        
        compind = conncomp(NTgraph);
        freq = hist(compind,1:max(compind));
        [~,largest] = max(freq);
        
        keepind = find(compind==largest);
        
        NT.keepNodes(keepind);
        
    end
    
    function [NTgraph,A,edgeN2G, edgeG2N] = makeGraph(NT)
        % convert network into a graph structure (via adjacency matrix)
        % edgeNT2G = for each network edge, corresponding edge index in
        % graph
        % edgeG2NT = for each graph edge, corresponding index in network

        % set up a connectivity matrix
        A = zeros(size(NT.nodepos,1),size(NT.nodepos,1));
        for ec = 1:size(NT.edgenodes,1)
            A(NT.edgenodes(ec,1),NT.edgenodes(ec,2)) = NT.edgelens(ec);
            A(NT.edgenodes(ec,2),NT.edgenodes(ec,1)) = NT.edgelens(ec);
        end
        % make into a graph structure
        NTgraph = graph(A);
        
        edgenodes = NT.edgenodes; % make earlier node first
        for ec = 1:NT.nedge
            edgenodes(ec,:) = sort(NT.edgenodes(ec,:));
        end

        [tmpN,Nsort] = sortrows(edgenodes,[1,2]);
        [tmpG,Gsort] = sortrows(NTgraph.Edges.EndNodes,[1,2]);

        edgeN2G(Nsort) = Gsort;
        edgeG2N(Gsort) = Nsort;
    end
    
    function [NTgraph,allcoords] = makeGraphEdgePath(NT)
        % create a graph structure including the points along the edge
        % paths
        
        %%
        % number of path points along each edge
        edgenpt = cellfun(@(x) size(x,1), NT.edgepath);
        % total number of points in new graph
        totpts = NT.nnode + sum(edgenpt-2);
        A = zeros(totpts,totpts);
        %%
        lastind = NT.nnode;
        
        allcoords = NT.nodepos;
        
        for ec = 1:NT.nedge
            n1 = NT.edgenodes(ec,1); n2 = NT.edgenodes(ec,2);
            
            newpts = NT.edgepath{ec}(2:end,:);
            % indices for new points
            newind = lastind+1:lastind+size(newpts,1);
            
            % keep track of all coordinates;
            allcoords(newind,:) = newpts;
            
            % connect to initial node
            A(newind(1),n1) = 1;
            A(n1,newind(1)) = 1;
            
            for pc = 2:size(newpts,1)
                A(newind(pc),newind(pc-1)) = 1;
                A(newind(pc-1),newind(pc)) = 1;
            end
            
            % connect to final node
            A(newind(end),n2) = 1;
            A(n2,newind(end)) = 1;
            
            lastind = lastind + size(newpts,1);
        end
        
        % make into a graph structure
        NTgraph = graph(A);
    end
    
    function edgepaths = nodepaths2edgepaths(NT,paths)
        % convert a series of node-to-node paths (eg: generated by
        % cyclebasis) to a series of directed edge paths 
        % edge index is negative if edge is traversed in reverse
              
        edgepaths = cell(size(paths));
        
        for p = 1:length(paths)
            path = paths{p};
            
            pathlen = length(path);
            path = [path path(1)];
            edgepath = zeros(1,pathlen);
            for cc = 1:pathlen
                nc1 = path(cc);
                nc2 = path(cc+1);
                ind = find(NT.nodenodes(nc1,:)==nc2);
                edge = NT.nodeedges(nc1,ind);
                % decide on edge direction
                if NT.edgenodes(edge,1)==nc1
                    edgepath(cc) = edge;
                elseif NT.edgenodes(edge,2)==nc1
                    edgepath(cc) = -edge;
                else
                    error('Problem with edge node connection')
                end
            end    
            
            edgepaths{p} = edgepath;
        end
    end
    
    function setLoops(NT)
        % calculate  independent loops in the network
        % stored in NT.loops as paths along edges
        % get adjacency matrix
        
        adjmat = zeros(NT.nnode,NT.nnode);
        for nc = 1:NT.nnode
            for cc = 1:NT.degrees(nc)
                nc2 = NT.nodenodes(nc,cc);
                adjmat(nc,nc2) = 1;
                adjmat(nc2,nc) = 1;
            end
        end
        
        %G = NT.network2Graph()               
        paths = cyclebasis(adjmat,'path');
        %
        edgepaths = NT.nodepaths2edgepaths(paths);
        
        NT.loops = edgepaths;
    end
    
    function setCumEdgeLen(NT, whichedges,setedgelens)
        %this fucntion will add a index to each pathcoord
        %  stores in NT.cumedgelen
        
        if (~exist('whichedges','var'))
            whichedges = 1:length(NT.edgepath);            
        end
        
        if (~exist('setedgelens','var'))
            setedgelens = 0;
        end
        
        for ec = whichedges
            path= NT.edgepath{ec};
            
            pathdiffs = path(2:end,:) - path(1:end-1,:);
            pathlens = sqrt(sum(pathdiffs.^2,2));
            
            NT.cumedgelen{ec} = [0,cumsum(pathlens')];
            
            if (setedgelens)
                % reset edge lens from cumulative               
                NT.edgelens(ec) = NT.cumedgelen{ec}(end);                
            end
        end
        
        if (setedgelens)
            % reset edge lens from cumulative
            NT.edgelens = NT.edgelens(1:NT.nedge);            
        end
    end
    
    function [G,A] = network2Graph(NT)
        % convert from our network object to a matlab graph object
        % returns graph G and adjacency matrix A
        
        % get adjacency matrix
        A = zeros(NT.nnode,NT.nnode);
        for ec = 1:NT.nedge
            nodes = NT.edgenodes(ec,:);
            A(nodes(1),nodes(2)) = 1;
            A(nodes(2),nodes(1)) = 1;
        end
        
        % specify edges and weights (edge lengths)
        G = graph(NT.edgenodes(:,1),NT.edgenodes(:,2),NT.edgelens);      
    end
    
     function interpolateEdgePaths(NT,npt,whichedges)
        % for a network object, get edge paths as a linear interpolation
        % between the end points of each edge
        % WARNING: this resets edge lengths to be the distances between end points
        % WARNING: should fix at some point to have similar step sizes in
        % interpolation (different npt for different length edges)
        
        % npt = number of points to interpolate along each edge.
        
        if NT.nedge==0
            % no edges to interpolate
            return
        end

        if (~exist('whichedges','var'))
            whichedges = 1:NT.nedge;
        end
        if (size(whichedges,1)>size(whichedges,2))
            whichedges = whichedges';
        end
        
        ipts = linspace(0,1,npt)'; % interpolation points
        for ec = whichedges
            n1 = NT.edgenodes(ec,1); n2 = NT.edgenodes(ec,2);
            p1 = NT.nodepos(n1,:); p2 = NT.nodepos(n2,:);
            
            NT.edgepath{ec} = p1 + ipts*(p2-p1);
            NT.edgelens(ec) = norm(p2-p1);
        end
        
     end
   
     function breakEdge(NT,ectarget,breakfrac,dosetup)
         % break up an edge in the network, creating a new degree-2 node along it
         % new node is located at fraction breakfrac along the edge
        % WARNING: does not currently save edgewidths

         if (~exist('dosetup','var'))
             dosetup = 1;
         end
         
        % get new node position
        shiftlen = NT.edgelens(ectarget)*breakfrac;
        newpos = interp1(NT.cumedgelen{ectarget},NT.edgepath{ectarget},shiftlen);
        % index that is right below the desired length
        tmp = find(NT.cumedgelen{ectarget}<shiftlen);        
        ind = max(tmp);
        %[~,ind] = min(abs(NT.cumedgelen{ectarget}-shiftlen));
        
        nnode = NT.nnode; nedge = NT.nedge;
        NT.nodepos(nnode+1,:) = newpos;        
        % add a new edge
        NT.edgenodes(end+1,:) = [NT.edgenodes(ectarget,2),nnode+1];
        NT.edgewidth{end+1} = [];
        % replace previous edge with one that goes to new node
        NT.edgenodes(ectarget,2) = nnode+1;
        
        if (~isempty(NT.nodevals))
            NT.nodevals(nnode+1) = 0;
        end
        if (~isempty(NT.nodelabels))
            NT.nodelabels{nnode+1} = '';
        end

        if (dosetup)
            % reset network
            NT.setupNetwork();

            % using old edgepath set new edgepaths as two pieces of original
            epath = NT.edgepath{ectarget};

            NT.edgepath{nedge+1} = flip([newpos; epath(ind+1:end,:)]);

            % in case newpos = a point along edgepath
            if (epath(ind,:)==newpos)
                NT.edgepath{ectarget}= epath(1:ind,:);
            else
                NT.edgepath{ectarget}= [epath(1:ind,:); newpos];
            end

            % set cumulative edge length
            NT.setCumEdgeLen();

            % then update edgelens (TODO should be able to do this in call to
            % cumedgelen....)
            NT.edgelens(nedge+1) = NT.cumedgelen{nedge+1}(end);
            NT.edgelens(ectarget) = NT.cumedgelen{ectarget}(end);

            if (~isempty(NT.edgevals))
                NT.edgevals(nedge+1) = {};
            end
            if (~isempty(NT.edgewidth))
                NT.edgewidth{nedge+1} = {};
            end

        end
     end
     
%      function breakEdge(NT,ectarget,breakfrac,dosetup)
%         % break up an edge in the network, creating a new degree-2 node 
%         % along it. new node is located at the edgepath point that is 
%         % closest to fraction breakfrac along edge
%         
%         if (~exist('dosetup','var'))
%             dosetup = 1;
%         end
% 
%         % get new node position
%         shiftlen = NT.edgelens(ectarget)*breakfrac;
%         [~,ind] = min(abs(NT.cumedgelen{ectarget}-shiftlen));
%         newpos = NT.edgepath{ectarget}(ind,:);
% 
%         nnode = NT.nnode; nedge = NT.nedge;
%         NT.nodepos(nnode+1,:) = newpos;
% 
%         % add a new edge
%         NT.edgenodes(end+1,:) = [nnode+1,NT.edgenodes(ectarget,2)];
%         % replace previous edge with one that goes to new node
%         NT.edgenodes(ectarget,2) = nnode+1;
% 
%         if (dosetup)
%             % reset network
%             NT.setupNetwork();
% 
%             % using old edgepath set new edgepaths as two pieces of original
%             epath = NT.edgepath{ectarget};
%             NT.edgepath{nedge+1} = epath(ind:end,:);
%             NT.edgepath{ectarget}= epath(1:ind,:);
% 
%             % set cumulative edge length
%             NT.setCumEdgeLen();
% 
%             % then update edgelens (TODO should be able to do this in call to
%             % cumedgelen....)
%             NT.edgelens(nedge+1) = NT.cumedgelen{nedge+1}(end);
%             NT.edgelens(ectarget) = NT.cumedgelen{ectarget}(end);
%         end
% 
%      end
     
         function outputPDB(NT,outfile,scl,sclval)
        %% output network as pdb formatted file
        % tracing along edge paths
                
        if (NT.dim~=3)
            error('pdb output only allowed for networks defined in 3D')
        end
        
        if (~exist('scl','var'))
            scl = 10;
        end
        if (~exist('sclval','var'))
            sclval = 1;
        end
        
        if (isempty(NT.edgepath))
            %NT.outputPDB_old(NT,outfile,scl);
            NT.interpolateEdgePaths(2);           
        end
        
        
        of = fopen(outfile,'w');
        
        % assign an index to each bead along the path
        % first in the index list are the nodes
        ct = NT.nnode;
        for ec = 1:NT.nedge
            edgepath = NT.edgepath{ec};
            beadind{ec} = [];
            for pc = 2:size(edgepath,1)-1
                ct = ct+1;
                beadlist(ct,:) = [ec pc];                
                beadind{ec}(pc) = ct;
            end
        end
               
        connections = zeros(ct,max(NT.degrees));
        ncon = zeros(ct,1);
        
        % output actual network nodes
        for nc = 1:NT.nnode           
            name = 'N';            
            
            if (~isempty(NT.edgevals))
                val = mean(NT.edgevals(NT.nodeedges(nc,1:NT.degrees(nc))));
            else
                val = 1;
            end
            
            
            pos = scl*NT.nodepos(nc,:);
            fprintf(of,'HETATM%5d%5s SSN X   0    %8.3f%8.3f%8.3f%6.2f%6.2f%13s\n',...
                nc,name,pos,1,val*sclval,'C');
            
            % keep track of all connections for this node
            connections(nc,1) = nc;            
            for ecc = 1:NT.degrees(nc)
                ec = NT.nodeedges(nc,ecc);
                if (isempty(beadind{ec}))
                    % connect to other nodes
                    if (NT.edgenodes(ec,1)==nc)
                        connections(nc,1+ecc) = NT.edgenodes(ec,2);
                    elseif (NT.edgenodes(ec,2)==nc)
                        connections(nc,1+ecc) =NT.edgenodes(ec,1);
                    else
                        error('bad network structure')
                    end
                else
                    % connect to nodes along the paths
                    if (NT.edgenodes(ec,1)==nc)
                        connections(nc,1+ecc) =beadind{ec}(2);
                    elseif(NT.edgenodes(ec,2)==nc)
                        connections(nc,1+ecc) =beadind{ec}(end);
                    else
                        error('bad network structure')
                    end
                end
                ncon(nc) = NT.degrees(nc);
            end
        end
        
        bct = NT.nnode;
        
        % output edge paths
        for ec = 1:NT.nedge
            name = 'EP';
            n1 = NT.edgenodes(ec,1); n2 = NT.edgenodes(ec,2);
            
            if (~isempty(NT.edgevals))
                val = NT.edgevals(ec);
            end
            
            edgepath = NT.edgepath{ec};
            for cc = 2:size(edgepath,1)-1  
                ind = beadind{ec}(cc);
                pos = scl*edgepath(cc,:);                               
                
                fprintf(of,'HETATM%5d%5s SSN X   0    %8.3f%8.3f%8.3f%6.2f%6.2f%13s\n',...
                         ind,name,pos,1,val*sclval,'C');
                
                if (cc==2)
                    con1 = n1;
                else
                    con1 = beadind{ec}(cc-1);
                end
                if (cc==size(edgepath,1)-1)
                    con2 = n2;
                else
                    con2 = beadind{ec}(cc+1);                    
                end
                
                connections(ind,1:3) = [ind con1 con2];
                ncon(ind) = 2;
            end            
        end
        
        % output connections
        for cc = 1:ct            
            constr = sprintf('%5d',connections(cc,1:ncon(cc)+1));
            fprintf(of,'CONECT%s\n',constr);
        end
        
        fclose(of);
     end

    
     
     function [nodeplotH,edgeplotH] = plotNetwork(NT,options)
         
         opt.labels = 0;
         opt.Parent = gca; % parent axes to plot on
         opt.nodeplotopt = {'filled'};
         opt.nodesize = 20;
         opt.edgestyle = '-';
         opt.edgecolor = [0 0 0];
         opt.nodecolor = [0 0 1];
         opt.plotnodes = 1:NT.nnode;
         opt.plotedges = 1;
         opt.edgeplotopt = {'LineWidth',.5};
         % plot curved paths instead of straight edges
         opt.plotedgepath = 1;
         % show data tips as edge or node index
         opt.datatipindex = false;
         % scaling factor
         opt.scl = 1;                  
         
         if (exist('options','var'))
             
             if (isfield(options,'plotoverimage') & options.plotoverimage)
                 % reset defaults for plotting over a BW image
                 opt.nodesize =20;
                 opt.nodecolor = [1 0 0];                 
                 opt.edgeplotopt = {'LineWidth',1,'Color','g'};
            end
             
             opt =copyStruct(options,opt,1);
         end
                      
         
         if (length(opt.nodesize)==1)
             opt.nodesize = opt.nodesize*ones(NT.nnode,1);
         end
         if (size(opt.nodecolor,1)==1)
             opt.nodecolor = opt.nodecolor.*ones(NT.nnode,3);
         end
         if (size(opt.edgecolor,1)==1)
             opt.edgecolor = opt.edgecolor.*ones(NT.nedge,3);
         end
         
         dim = NT.dim;
         nodepos = NT.nodepos; 
         edgenodes = NT.edgenodes;
         scl = opt.scl;
         
         % if (size(opt.nodecolor,1)==1)
         %     opt.nodecolor = repmat(opt.nodecolor,size(nodepos,1),1);
         % end
         
         dttemplateset = false;
         if (opt.plotedges)
             for ec = 1:size(edgenodes,1)
                 if (~isempty(NT.edgepath) & opt.plotedgepath)
                     % plot curved paths of the edges
                     path = NT.edgepath{ec};
                     if (isempty(path))
                         n1 = NT.edgenodes(ec,1); n2 = NT.edgenodes(ec,2);
                         path = [NT.nodepos(n1,:); NT.nodepos(n2,:)];
                     end
                     
                     if (dim==2)
                         edgeplotH(ec) = plot(path(:,1)*scl,path(:,2)*scl,opt.edgestyle,'Color',opt.edgecolor(ec,:),opt.edgeplotopt{:},'Parent',opt.Parent);
                     else
                         edgeplotH(ec) = plot3(path(:,1)*scl,path(:,2)*scl,path(:,3)*scl,opt.edgestyle,'Color',opt.edgecolor(ec,:),opt.edgeplotopt{:});
                     end
                     
                     % label the edge graphics object with the
                     % corresponding index
                     % only label if the edgepath is > 2 points long (consists not
                     % just of endpoints)
                     if (~isempty(NT.edgepath) & opt.plotedgepath & opt.datatipindex)
                          edgeplotH(ec).addprop('edgeind');
                             edgeplotH(ec).edgeind = ec;
                         if (size(NT.edgepath{ec},1)>2)                                                         
                             dttemplate = edgeplotH(ec).DataTipTemplate;
                             dttemplate.FontSize=6;
                             dttemplate.DataTipRows(1).Value = ec*ones(size(NT.edgepath{ec},1),1);
                             dttemplate.DataTipRows(1).Label = '';
                             dttemplate.DataTipRows(2:end) = [];
                             dttemplateset = true;
                         else
                             edgeplotH(ec).PickableParts = 'none';
                         end
                     end
                     
                     % turn off datatips for edge paths
                     %edgeplotH(ec).PickableParts = 'none';
                 else
                     p1 = edgenodes(ec,1); p2 = edgenodes(ec,2);
                     if (dim==2)
                         plot(nodepos([p1,p2],1)*scl,nodepos([p1,p2],2)*scl,'Color',opt.edgecolor(ec,:),opt.edgeplotopt{:},'Parent',opt.Parent);
                     else
                         plot3(nodepos([p1,p2],1)*scl,nodepos([p1,p2],2)*scl,nodepos([p1,p2],3)*scl,'Color',opt.edgecolor(ec,:),opt.edgeplotopt{:});
                     end
                 end
                 hold all
             end
             axis equal
         end
         if (~isempty(opt.plotnodes))
             if (dim==2)
                 nodeplotH = scatter(nodepos(opt.plotnodes,1)*scl,nodepos(opt.plotnodes,2)*scl,opt.nodesize(opt.plotnodes),...
                     opt.nodecolor(opt.plotnodes,:),opt.nodeplotopt{:},'Parent',opt.Parent);
             else
                 nodeplotH = scatter3(nodepos(opt.plotnodes,1)*scl,nodepos(opt.plotnodes,2)*scl,nodepos(opt.plotnodes,3)*scl,...
                     opt.nodesize(opt.plotnodes),opt.nodecolor(opt.plotnodes,:),opt.nodeplotopt{:});
             end
             axis equal
             hold all
         end
         
         % set data tip properties
         if (opt.datatipindex & ~isempty(opt.plotnodes))
             dt = nodeplotH.DataTipTemplate;
             dt.DataTipRows(1).Value = 1:NT.nnode;
             dt.DataTipRows(1).Label = '';
             dt.DataTipRows(2:end) = [];
             dt.FontSize=6;
         end
         
         
         if (opt.labels)
             for pc = 1:length(nodepos)
                 text(nodepos(pc,1),nodepos(pc,2),sprintf('%d',pc))
             end
         end
         hold off
     end
     
     function interplen = reinterpolateEdgePaths(NT,dxwant,options)
         % reinterpolate the edge paths to have points
         % at *approximately* the desired spacing
         
         opt = struct();
         opt.interpmethod = 'linear';
         
         if (exist('options','var'))
             opt = copyStruct(options,opt);
         end
         
         if (isempty(NT.edgepath))
             error('no edge paths computed')
         end
         
         NT.setCumEdgeLen();                  
         interplen = {};
         for ec = 1:NT.nedge
             path = NT.edgepath{ec};
             cumlen = NT.cumedgelen{ec};
             
             % interpolate to this many equal-sized segments (roughly)
             % will have slightly smaller segments if not integer multiple
             nseg = ceil(cumlen(end)/dxwant);
             interpx = linspace(0,cumlen(end),nseg+1);
             interplen{ec} = interpx;
             
             newpts = interp1(NT.cumedgelen{ec},path,interpx,opt.interpmethod);
             NT.edgepath{ec} = newpts;
         end
         
         % recalculate cumulative lengths and total edge lengths
         NT.setCumEdgeLen(1:NT.nedge,true)
     end
     
     function setEdgeLens(NT,euclidean)
         % recalculate edge lengths, based on euclidean distance or
         % cumedgelens
         if (~exist('euclidean'))
             euclidean = false;
         end
         
         NT.edgelens = zeros(NT.nedge,1);
         
         if (euclidean)
             for ec = 1:NT.nedge
                 n1 = NT.edgenodes(ec,1); n2 = NT.edgenodes(ec,2);
                 NT.edgelens(ec) = norm(NT.nodepos(n2,:)-NT.nodepos(n1,:));
             end
         else
             if (isempty(NT.cumedgelen))
                 NT.setCumEdgeLen(1:NT.nedge,true);
             else
                 for ec = 1:NT.nedge
                     NT.edgelens(ec) = NT.cumedgelen{ec}(end);
                 end
             end
         end
     end
     
     function [xy,edgeind,edgepos,distance] = getNearestEdgePoint(NT,pts)
         % for a set of points (pts, size nxdim) find the nearest point
         % along the edge paths of the network
         % returns: 
         % xy = coordinates of the nearest points along interpolated edge paths
         % edgeind = indices of nearest edge to each point
         % edgepos = interpolated fractional distance along each edge to the point
         % distance = distance from original point to the mapped one
         % uses linear interpolation
         
         if (isempty(NT.edgepath))
             error('need to set edge path first')
         end
         
         distance = inf*ones(size(pts,1),1);
         for ec = 1:NT.nedge
             curvexy = NT.edgepath{ec};
             [xytmp,dist,t_a] = distance2curve(curvexy,pts);
             
             % which indices to update
             ind = find(dist<distance);
             
             if (~isempty(ind))
                 xy(ind,:) = xytmp(ind,:);
                 distance(ind) = dist(ind);
                 edgepos(ind) = t_a(ind);
                 edgeind(ind) = ec;
             end
         end
         
     end
     
     function fixEdgePathOrientation(NT)
         % reorient all edge paths to correctly match up to edge directions
         % also recalculates cumulative edge lengths
         for ec = 1:NT.nedge
             n1 = NT.edgenodes(ec,1);
             d1 = norm(NT.edgepath{ec}(1,:)-NT.nodepos(n1,:));
             d2 = norm(NT.edgepath{ec}(end,:)-NT.nodepos(n1,:));
             
             if (d1>d2)
                 NT.edgepath{ec} = flipud(NT.edgepath{ec});
             end
         end
         
         NT.setCumEdgeLen()
         
     end
     
     function mergednodes = mergeAllEdgePaths(NT,outedgeonly)
         %% remove all degree 2 nodes from the network, converting to longer
         % edges with defined edgepath instead
         % if outedgeonly is set (default=false), then only trace out
         % outgoing edges for each node. 
         % this ensures directed trees stay directed and edge widths
         % for short edges (no intermediate nodes) are set by downstream
         % end-node width
         % mergednodes stores all original nodes that went into the edge in
         % the final network
         % WARNING: will always keep root node, even if it has degree 2

        
         if (~exist('outedgeonly','var'))
             outedgeonly = false;
         end
         
         % find all nodes not of degree 2
         degnot2 = find(NT.degrees~=2);
         
         dokeep = true(1,NT.nnode);
         havechecked = false(1,NT.nedge);
         
         % which nodes got merged away for each edge
         mergednodes = {};
         for ec = 1:NT.nedge
             mergednodes{ec} = [];
         end
         
         % start with root node if set
         if (~isnan(NT.rootnode))
             startnodes = [NT.rootnode;degnot2(degnot2~=NT.rootnode)];
         else
             startnodes = degnot2;
         end
         
         for ncc = 1:length(startnodes)
             nc = startnodes(ncc);
                         
             for bc =1:NT.degrees(nc)
                 nclist= [nc]; edgelist = []; edgedir = [];
                 
                 % start new edge path
                 path = NT.nodepos(nc,:);

                 % go to each neighbor of this node
                 ncnext = NT.nodenodes(nc,bc);
                 
                 ec = NT.nodeedges(nc,bc);
                 % already got rid of nodes in this direction?
                 if (havechecked(ec)); continue; end                    
                 
                 nclist(end+1) = ncnext;                 
                 degnext = NT.degrees(ncnext);
                 havechecked(ec) = true;
                 
                 % list of original edges merging into the big one
                 ecnext = NT.nodeedges(nc,bc);
                 edgelist(end+1) = ecnext;                 
                 % direction of original edges
                 if (NT.edgenodes(ecnext,1) ==nc)
                     edgedir(end+1) = 1; % outgoing edge

                     if (~isempty(NT.edgepath))
                         newpath = NT.edgepath{ecnext}(2:end,:);
                         path = [path; newpath];
                     end
                 else
                     edgedir(end+1) = -1; % incoming edge


                     if (~isempty(NT.edgepath))
                         newpath = NT.edgepath{ecnext}(1:end-1,:);
                         path = [path; flipud(newpath)];
                     end
                 end

                 if (outedgeonly & edgedir(end)==-1)
                     % only following along outgoing edges
                     havechecked(ec) = false;
                     continue
                 end
                 
                 if (degnext~=2)
                     % have reached a junction or terminus
                     continue
                 end
                 
                  % keep building up edge
                 while degnext==2
                     % keep building up edge                    
                     bc2 = find(NT.nodenodes(ncnext,1:degnext)~=nclist(end-1));
                     if (length(bc2)>1)
                         error('bad connection')
                     end
                     % list of original edges merging into the big one
                     ecnext = NT.nodeedges(ncnext,bc2);
                     edgelist(end+1) = ecnext;
                     havechecked(ecnext) = true;
                     
                     ncnext = NT.nodenodes(ncnext,bc2);
                     degnext = NT.degrees(ncnext);
                     nclist(end+1) = ncnext; 
                    
                     % direction of original edges
                     if (NT.edgenodes(ecnext,1)==ncnext)
                         edgedir(end+1) = -1; % reverse direction

                         if (~isempty(NT.edgepath))
                             newpath = NT.edgepath{ecnext}(1:end-1,:);
                             path = [path; flipud(newpath)];
                         end
                     else
                         edgedir(end+1) = 1;

                         if (~isempty(NT.edgepath))
                             newpath = NT.edgepath{ecnext}(2:end,:);
                             path = [path; newpath];
                         end
                     end

                
                 end
                 
                 % make a path for this edge
                % path = NT.nodepos(nclist,:);
                 
                 % make new long edge
                 ecnew = size(NT.edgenodes,1)+1;
                 NT.edgenodes(ecnew,:) = [nclist(1) nclist(end)];
                 
                 % make a path for this edge
                 %NT.edgepath{ecnew} = NT.nodepos(nclist,:);
                 if (~isempty(NT.edgepath))
                     NT.edgepath{ecnew} = path;
                 end
                 
                 NT.setCumEdgeLen(ecnew);
                 NT.edgelens(ecnew) = NT.cumedgelen{ecnew}(end);
                 if (~isempty(NT.edgevals))
                     if (iscell(NT.edgevals))
                         NT.edgevals{ecnew} = [];
                     else
                        NT.edgevals(ecnew) = 0;
                     end
                 end
                 
                 % copy over width measurements
                 if (~isempty(NT.edgewidth))
                     NT.edgewidth{ecnew} = [];
                     for ecc = 1:length(edgelist)
                         ec = edgelist(ecc);
                         
                         start = NT.cumedgelen{ecnew}(ecc);
                         if (edgedir(ecc)>0)
                             addwidths =  NT.edgewidth{ec};
                             if (~isempty(addwidths))
                                addwidths(:,2) = addwidths(:,2)+start;
                             end
                         else
                             addwidths =  flipud(NT.edgewidth{ec});
                             if (~isempty(addwidths))
                                addwidths(:,2) = NT.edgelens(ec)-addwidths(:,2); % flip edge direction
                                addwidths(:,2) = addwidths(:,2)+start;
                             end
                         end
                         NT.edgewidth{ecnew} = [ NT.edgewidth{ecnew};addwidths];
                     end
                 end
                 
                 % keep track of nodes that got merged away to make this
                 % new edge
                 mergednodes{size(NT.edgenodes,1)} = nclist;
                 
                 % mark intermediate nodes for removal
                 dokeep(nclist(2:end-1)) = false;
             end             
         end
         
         % update connectivity info
         NT.setupNetwork();
         
         % get rid of excess nodes
         keepnodes = find(dokeep);
         [mapold2new,mapnew2oldedge] = NT.keepNodes(keepnodes);
         
         mergednodes = mergednodes(mapnew2oldedge);
         
     end
     
     function scaleCoords(NT,scl)
        % multiply all spatial coordinates by scl
         disp('Scaling all coords, including edge widths')
        NT.nodepos = NT.nodepos*scl;
        for ec = 1:NT.nedge
            if (~isempty(NT.edgepath))
                NT.edgepath{ec} = NT.edgepath{ec}*scl;
            end
            if (~isempty(NT.cumedgelen))
                NT.cumedgelen{ec} =NT.cumedgelen{ec}*scl;
            end
        end
        NT.edgelens = NT.edgelens*scl;

        % update edge widths if they exist
         if (~isempty(NT.edgewidth))
             for ec = 1:NT.nedge
                NT.edgewidth{ec} = NT.edgewidth{ec}*scl; 
             end
         end
     end

     function scaleCoords3D(NT,scl)
         % scale all spatial coordinates by different amounts in the x, y,
         % z dimensions. 
         % scl(i) = um per pixel in the i dimension for i = 1,2,3
         
         for i = 1:size(NT.nodepos,2)
             NT.nodepos(:,i) = NT.nodepos(:,i)*scl(i);             
         end

         if (~isempty(NT.edgepath))
             for ec = 1:NT.nedge
                 if (~isempty(NT.edgepath{ec}))
                     for i = 1:size(NT.nodepos,2)
                         NT.edgepath{ec}(:,i) = NT.edgepath{ec}(:,i)*scl(i);
                     end
                 end
             end
         end

         % update edge widths if they exist
         if (~isempty(NT.edgewidth))
             if (abs(scl(1)-scl(2))>10*eps)
                 error('Currently cannot rescale edge widths if x and y scalings are different: %f', scl)
             end
             for ec = 1:NT.nedge
                NT.edgewidth{ec} = NT.edgewidth{ec}*scl(1); 
             end
         end

         % update all edge lengths
         NT.setCumEdgeLen(1:NT.nedge, true);
     end
end

end