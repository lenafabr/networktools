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
    edgelens % edge lengths
    nodelabels % string label for each node (permeability or reservoir)
    nodevals
    edgevals
    loops
    
    edgepath % coordinates of path along each edge
    cumedgelen % cumulative lengths along each edge
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
        
        edgevals
        
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
        end
            
    end
    
    function removeDoubleEdges(NT)
        %% get rid of all duplicate edges, so you only have one edge
        % between each pair of nodes
        badind = [];
        newedgenodes = [];
        edgenodes= NT.edgenodes
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
    
    function [mapold2new] = keepNodes(NT,keepind)        
        % truncate network to only keep the node indices given by keepind
        % remove any edges not between nodes in the new network.
        % keep only largest connected component
        % mapnew2old(i) = index of new node i in the old node list
        % mapold2new(i) = index of old node i in new node list
        
        [newnodepos,newedgenodes,mapold2new] = truncateNetworkNodes(keepind,NT.nodepos,NT.edgenodes);
        
        NT.nodepos = newnodepos;
        NT.edgenodes = newedgenodes;
        
        NT.setupNetwork()
                
        if (~isempty(NT.nodevals)); NT.nodevals = NT.nodevals(keepind); end
        if (~isempty(NT.edgevals)); NT.edgevals = NT.edgevals(keepind); end
     
    end
    
    function addNodes(NT,nodepos,connect)
        % add nodes to the network
        % connecting each to predefined old node indices
        nnew = size(nodepos,1);
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
    end
    
    function outputNetwork(NT,outfile,edgepaths,nodelabels)
        % output network structure to file
        % edgepaths: lists of cycles, in terms of edges to include as LOOP
        % structures
        % reservoirnodes: which nodes belong to a certain reservoir
        
        of = fopen(outfile,'w')
        fprintf(of,'%s\n','# file defining network structure')
        fprintf(of,'%s\n\n','# made with matlab NetworkObj output')
        
        fprintf(of,'%s \n','# list of node indices and xy positions and values')        
        nodefmtstring = ['NODE %d ' repmat(['%20.10f '],1,NT.dim+1) '\n'];
        nodelblfmtstring = ['NODE %d ' repmat(['%20.10f '],1,NT.dim+1) '%s \n'];
        %
        for pc = 1:size(NT.nodepos)
            
            if (isempty(NT.nodevals))
                val = 0.0;              
            else
                val = NT.nodevals(pc);                
            end
            
            if (exist('nodelabels','var'))
                fprintf(of,nodelblfmtstring,pc, NT.nodepos(pc,:), val,nodelabels{pc});
            elseif (~isempty(NT.nodelabels))
                fprintf(of,nodelblfmtstring,pc, NT.nodepos(pc,:), val,NT.nodelabels{pc});
            else
                fprintf(of,nodefmtstring,pc, NT.nodepos(pc,:), val);
            end
        end
        % edge information
        edgefmtstring = 'EDGE %d %d %d %20.10f\n';
        for ec = 1:size(NT.edgenodes)
            fprintf(of,edgefmtstring,ec, NT.edgenodes(ec,:), NT.edgelens(ec));
        end
        
        % independent cycles as a list of edges
        if (exist('edgepaths','var'))
            for pc = 1:length(edgepaths)
                fprintf(of,['LOOP ' sprintf('%d ',edgepaths{pc}) '\n'])
            end   
        elseif (~isempty(NT.loops))
            for pc = 1:length(NT.loops)
                fprintf(of,['LOOP ' sprintf('%d ',NT.loops{pc}) '\n'])
            end   
        end
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
    
    function [NTgraph,A] = makeGraph(NT)
        % convert network into a graph structure (via adjacency matrix)
        
        % set up a connectivity matrix
        A = zeros(size(NT.nodepos,1),size(NT.nodepos,1));
        for ec = 1:size(NT.edgenodes,1)
            A(NT.edgenodes(ec,1),NT.edgenodes(ec,2)) = 1;
            A(NT.edgenodes(ec,2),NT.edgenodes(ec,1)) = 1;
        end
        % make into a graph structure
        NTgraph = graph(A);
        
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
        
        numedge=length(NT.edgepath);%get cell size
        for num =whichedges%go to individual cell
            [seg,~]=size(NT.edgepath{num});%get path coord size
            dx(1)=0;%set a distance vector for x
            dy(1)=0;%set a distance vector for y
            dx(2:seg)=NT.edgepath{num}(2:seg,1)-NT.edgepath{num}(1:seg-1,1);
            dy(2:seg)=NT.edgepath{num}(2:seg,2)-NT.edgepath{num}(1:seg-1,2);
            if (NT.dim==3)
                dz(1) = 0;
                dz(2:seg) = NT.edgepath{num}(2:seg,3)-NT.edgepath{num}(1:seg-1,3);
                dl = sqrt(dx.^2+dy.^2+dz.^2);
            else
                dl = sqrt(dx.^2+dy.^2);
            end
            for segnum = 1:seg
                L= sum(dl(1:segnum));
                NT.cumedgelen{num}(segnum)=L;
            end
        end
        
        
        if (setedgelens)
            % reset edge lens from cumulative
            for ec = 1:NT.nedge
                NT.edgelens(ec) = NT.cumedgelen{ec}(end);
            end
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
        
        G = graph(A);                
    end
    
     function interpolateEdgePaths(NT,npt,whichedges)
        % for a network object, get edge paths as a linear interpolation
        % between the end points of each edge
        % WARNING: this resets edge lengths to be the distances between end points
        % WARNING: should fix at some point to have similar step sizes in
        % interpolation (different npt for different length edges)
        
        % npt = number of points to interpolate along each edge.
        
        if (~exist('whichedges','var'))
            whichedges = 1:NT.nedge;
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
        
         if (~exist('dosetup','var'))
             dosetup = 1;
         end
         
        % get new node position
        shiftlen = NT.edgelens(ectarget)*breakfrac;
        newpos = interp1(NT.cumedgelen{ectarget},NT.edgepath{ectarget},shiftlen);
        
        nnode = NT.nnode; nedge = NT.nedge;
        NT.nodepos(nnode+1,:) = newpos;        
        % add a new edge
        NT.edgenodes(end+1,:) = [NT.edgenodes(ectarget,2),nnode+1];
        % replace previous edge with one that goes to new node
        NT.edgenodes(ectarget,2) = nnode+1;
        
        if (dosetup)
        % reset network
        NT.setupNetwork();
        
        % interpolate paths on these new edges
        NT.interpolateEdgePaths(2,[ectarget,nedge+1]);
        NT.setCumEdgeLen(1);
        end
     end
     
     function outputPDB(NT,outfile,scl)
        %% output network as pdb formatted file
        % tracing along edge paths
                
        if (~exist('scl','var'))
            scl = 10;
        end
        
        if (isempty(NT.edgepath))
            NT.outputPDB_old(NT,outfile,scl);
            return 
        end
        
        
        of = fopen(outfile,'w');
        
        % assign an index to each bead along the path
        % first in the index list are the nodes
        ct = NT.nnode;
        for ec = 1:NT.nedge
            edgepath = NT.edgepath{ec};
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
            
            pos = scl*NT.nodepos(nc,:);
            fprintf(of,'HETATM%5d%5s SSN X   0    %8.3f%8.3f%8.3f%6.2f%6.2f%13s\n',...
                nc,name,pos,1,1,'C');
            
            % keep track of all connections for this node
            connections(nc,1) = nc;            
            for ecc = 1:NT.degrees(nc)
                ec = NT.nodeedges(nc,ecc);
                if (NT.edgenodes(ec,1)==nc)
                    connections(nc,1+ecc) =beadind{ec}(2);
                elseif(NT.edgenodes(ec,2)==nc)
                    connections(nc,1+ecc) =beadind{ec}(end);
                else
                    error('bad network structure')
                end
                ncon(nc) = NT.degrees(nc);
            end
        end
        
        bct = NT.nnode;
        
        % output edge paths
        for ec = 1:NT.nedge
            name = 'EP';
            n1 = NT.edgenodes(ec,1); n2 = NT.edgenodes(ec,2);
            
            edgepath = NT.edgepath{ec};
            for cc = 2:length(edgepath)-1  
                ind = beadind{ec}(cc);
                pos = scl*edgepath(cc,:);
                
                fprintf(of,'HETATM%5d%5s SSN X   0    %8.3f%8.3f%8.3f%6.2f%6.2f%13s\n',...
                         ind,name,pos,1,1,'C');
                
                if (cc==2)
                    con1 = n1;
                else
                    con1 = beadind{ec}(cc-1);
                end
                if (cc==length(edgepath)-1)
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
            constr = sprintf('%5d',connections(cc,1:ncon(cc)));
            fprintf(of,'CONECT%s\n',constr);
        end
        
        fclose(of);
     end
     
     function plotNetwork(NT,options)
         opt.labels = 0;
         opt.nodeplotopt = {'b','filled'};
         opt.nodesize = 20;
         %opt.nodecolor = [0 0 0];
         opt.plotnodes = 1:NT.nnode;
         opt.plotedges = 1;
         opt.edgeplotopt = {'MarkerSize',1};
         % plot curved paths instead of straight edges
         opt.plotedgepath = 1;
         
         if (exist('options','var'))
             opt =copyStruct(options,opt,1);
         end
         
         nodepos = NT.nodepos;
         edgenodes = NT.edgenodes;
         
         
         if (length(opt.nodesize)==1)
             opt.nodesize = opt.nodesize*ones(size(nodepos,1),1);
         end
                           
         dim = size(nodepos,2);
         % if (size(opt.nodecolor,1)==1)
         %     opt.nodecolor = repmat(opt.nodecolor,size(nodepos,1),1);
         % end
         
         if (opt.plotedges)
             for ec = 1:length(edgenodes)
                 if (opt.plotedgepath)
                     % plot curved paths of the edges
                     path = NT.edgepath{ec};
                     if (dim==2)
                         plot(path(:,1),path(:,2),'k',opt.edgeplotopt{:})
                     else
                         plot3(path(:,1),path(:,2),path(:,3),'k.-',opt.edgeplotopt{:})
                     end
                 else
                     p1 = edgenodes(ec,1); p2 = edgenodes(ec,2);
                     if (dim==2)
                         plot(nodepos([p1,p2],1),nodepos([p1,p2],2),'k',opt.edgeplotopt{:})
                     else
                         plot3(nodepos([p1,p2],1),nodepos([p1,p2],2),nodepos([p1,p2],3),'k',opt.edgeplotopt{:})
                     end
                 end
                 hold all
             end
             axis equal
         end
         if (~isempty(opt.plotnodes))
             %plot(nodepos(:,1),nodepos(:,2),'b.',opt.nodeplotopt{:})
             %scatter(nodepos(opt.plotnodes,1),nodepos(opt.plotnodes,2),opt.nodesize(opt.plotnodes),opt.nodecolor(opt.plotnodes)',opt.nodeplotopt{:})
             if (dim==2)
                 scatter(nodepos(opt.plotnodes,1),nodepos(opt.plotnodes,2),opt.nodesize(opt.plotnodes),opt.nodeplotopt{:})
             else
                 scatter3(nodepos(opt.plotnodes,1),nodepos(opt.plotnodes,2),nodepos(opt.plotnodes,3),...
                     opt.nodesize(opt.plotnodes),opt.nodeplotopt{:})
             end
             axis equal
             hold all
         end
         
         if (opt.labels)
             for pc = 1:length(nodepos)
                 text(nodepos(pc,1),nodepos(pc,2),sprintf('%d',pc))
             end
         end
         hold off
     end
     
     function reinterpolateEdgePaths(NT,dxwant)
         % reinterpolate the edge paths to have points
         % at *approximately* the desired spacing
         
         if (isempty(NT.edgepath))
             error('no edge paths computed')
         end
         
         NT.setCumEdgeLen();                  
         
         for ec = 1:NT.nedge
             path = NT.edgepath{ec};
             cumlen = NT.cumedgelen{ec};
             
             % interpolate to this many equal-sized segments (roughly)
             % will have slightly smaller segments if not integer multiple
             nseg = ceil(cumlen(end)/dxwant);
             interpx = linspace(0,cumlen(end),nseg+1);
             
             newpts = interp1(NT.cumedgelen{ec},path,interpx);
             NT.edgepath{ec} = newpts;
         end
         
         % recalculate cumulative lengths and total edge lengths
         NT.setCumEdgeLen(1:NT.nedge,true)
     end
end
end