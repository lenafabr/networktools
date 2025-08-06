% class definition for one of our network objects
% this mostly relies on Matlab's builtin digraph object
% with some additional fields and methods included

classdef NetworkGraphObj < matlab.mixin.Copyable   
    properties
        graph
        dim % spatial dimension
        % last unique index used for naming nodes
        lastnode=0 
    end
    
    methods
        function NT = NetworkGraphObj(fname,options)
            NT.graph = digraph;
            NT.dim = 0;
            
            if (exist('fname','var'))
                % read in network data from file            
                error('reading from file not set upyet')
            end               
        end
        
        function setupNetwork(NT,edgenodes,nodepos,edgepath)
            % set up network info from arrays of edge connectivites,
            % node positions, and (optionally) paths for the edges
            % edgenodes = m x 2 array ofconnected pairs of nodes
            % nodepos = n x dim array of node positions
            % edgepath = cell array with p x dim arrays of edge paths                       
            
            names = {};
            for nc = 1:size(nodepos,1)
                NT.lastnode = NT.lastnode + 1;
                names{nc} = sprintf('n%d',NT.lastnode);
            end            
            
            NT.graph = digraph(edgenodes(:,1),edgenodes(:,2),[],names);
            
            NT.dim = size(nodepos,2);
            
            NT.graph.Nodes.pos = nodepos;            
            
            if (exist('edgepath','var'))
                NT.graph.Edges.edgepath = edgepath;
            end
            
            NT.setCumEdgeLen(1:NT.graph.numedges,true);
                       
        end
            
        function [graphplotH,edgeplotH] = plotNetwork(NT,options)
            G = NT.graph;
            
            % default options
            opt.edgeplotopt = {'LineWidth',1};
            % show straight edge directions
            opt.showEdgeDir = 0;
            % edge color (all black by default; can be preset for each
            % edge)
            opt.edgecolor = [0 0 0];
            % parent axes on which to plot;
            opt.Parent = [];           
            
            if (size(opt.edgecolor,1)==1)
                opt.edgecolor = opt.edgecolor.*ones(G.numedges,3);
            end
            % color of nodes
            opt.nodecolor = [0 0 1];
            % size ofnodes
            opt.nodesize = 5;            
            opt.edgelabels = [];
         
            if (exist('options','var'))
             opt =copyStruct(options,opt,1);
            end
            
            if (isempty(opt.Parent)); opt.Parent = gca; end
            
            if (NT.dim>2)
                error('Plotting currently only set up for 2D networks')
            end
                        
            
            for ec = 1:G.numedges
                path = G.Edges.edgepath{ec};
                edgeplotH(ec) = plot(path(:,1),path(:,2),'Color',opt.edgecolor(ec,:),opt.edgeplotopt{:},...
                    'Parent',opt.Parent);
                edgeplotH(ec).PickableParts = 'none';
                hold all
            end            
            graphplotH = G.plot('XData',G.Nodes.pos(:,1),'YData',G.Nodes.pos(:,2),...,
                'EdgeAlpha',opt.showEdgeDir,...
                'NodeColor',opt.nodecolor,'MarkerSize',opt.nodesize,'Parent',opt.Parent);
            graphplotH.DataTipTemplate.FontSize = 6;
            
            % label edge numbers
            if (~isempty(opt.edgelabels))
                labeledge(graphplotH,1:G.numedges,opt.edgelabels)
            end
            
            % make data tips only show node number, not extra info
            graphplotH.DataTipTemplate.DataTipRows(2:end) = [];
            hold off
        end
        
        function setCumEdgeLen(NT, whichedges,setedgelens)
        % calculate cumulative lengths along the edge paths
        % whichedges: optionally decide which edges to calculate
        % setedgelens: optionally, set the edgelen field as well
        
        G = NT.graph; Edges = G.Edges;
        if (~ismember('edgepath',Edges.Properties.VariableNames))
            error('edge paths not set')
        end
        
        if (~exist('whichedges','var'))
            whichedges = 1:G.numedges;            
        end
        
        if (~exist('setedgelens','var'))
            setedgelens = false;
        end
        
        
        for ec = whichedges
            path= Edges.edgepath{ec};
            
            pathdiffs = path(2:end,:) - path(1:end-1,:);
            pathlens = sqrt(sum(pathdiffs.^2,2));
            
            G.Edges.cumedgelen{ec} = [0,cumsum(pathlens')];
            
            if (setedgelens)
                % reset edge lens from cumulative               
                edgelens(ec) = G.Edges.cumedgelen{ec}(end);                
            end
        end
               
        if (setedgelens)
            % reset edge lens from cumulative
            G.Edges.edgelens(whichedges) = edgelens(whichedges);            
        end
        
        NT.graph = G;
    end
        
    function setEdgeLens(NT,euclidean)
         % recalculate edge lengths, based on euclidean distance or
         % cumedgelens
         if (~exist('euclidean'))
             euclidean = false;
         end
         
         G = NT.graph; Edges = G.Edges; Nodes = G.Nodes;
         
         edgelens = zeros(G.numedges,1);
         
         if (euclidean)
             for ec = 1:G.numedges
                 n1 = Edges.EndNodes(ec,1); n2 = Edges.EndNodes(ec,2);                 
                 edgelens(ec) = norm(Nodes.pos(n2,:)-Nodes.pos(n1,:));
             end
         else
             if (~ismember('cumedgelen',Edges.Properties.VariableNames))
                 % calculate cumulate lengths if not already set
                 NT.setCumEdgeLen(1:NT.nedge,true);
             else % set edge length based on cumulative lens
                 for ec = 1:G.numedges
                     edgelens(ec) = Edges.cumedgelen{ec}(end);
                 end
             end
         end
         
         NT.graph.Edges.edgelens = edgelens;
    end
        
     function removeDoubleEdges(NT)
        %% get rid of all duplicate edges, so you only have one edge
        % between each pair of nodes
        
        NT.graph = simplify(NT.graph);      
     end
    
     function keepNodes(NT,keepind)
         % truncate network to keep only the node indices given by keepind
         % this is only here for backward compatibility with NetworkObj
         NT.graph = subgraph(NT.graph,keepind);
     end
     
     function keepEdges(NT,edgeind)
         % keep only the given edges in a graph
         % this is only here for backward compatibility with NetworkObj
         willkeepedge = false(1,NT.graph.numedges);
         willkeepedge(edgeind) = true;
         
         NT.graph = rmedge(NT.graph,find(~willkeepedge));
     end
     
     function rmEdges(NT,edgeind)
         % remove specific edges
         NT.graph = rmedge(NT.graph,edgeind);
     end
     
     function newind = addEdges(NT,n1,n2,paths)
         % add edge between two nodes
         % n1 and n2 are vectors of node pairs
         % with the specified paths (paths is a cell array of path coords)
         % if no paths given, just do a direct line between the 2 nodes
         
         NT.graph = addedge(NT.graph,n1,n2);
         newind = findedge(NT.graph,n1,n2)';     
         
         if (exist('paths','var'))
            NT.graph.Edges.edgepath(newind) = paths;
         else
             NT.interpolateEdgePaths(2,newind);
         end
         
         % set length and cumulative length
         NT.setCumEdgeLen(newind,true);
     end
     
     function newedgeind = addNodes(NT,nodepos,connect,edgepaths)
         % add nodes to the network graph
         % connect each one to a set of other node indices (incoming edges
         % to the new node)
         % optionally, preset the edge paths  (otherwise, interpolate) 
         % edgepath is a cell array listing the edge paths for each new
         % edge, for each of the new nodes in order
         
         nnode = NT.graph.numnodes;
         nedge = NT.graph.numedges;
         
         nnew = size(nodepos,1);
         %newind = nnode+1:nnode+nnew; 
              
         newnames = cell(nnew,1);
         for nc = 1:nnew
             NT.lastnode = NT.lastnode+1;
             newnames{nc} = sprintf('n%d',NT.lastnode);
         end
         
         % add connections for the nodes
         newedges = {};
         if (exist('connect','var'))
             for pc = 1:nnew
                 nedge = length(connect{pc});
                 for ec= 1:nedge
                     newedges(end+1,:) = {connect{pc}{ec} newnames{pc}};
                 end
             end
             nnewedge = size(newedges,1);
         end                 
         
         nodetab = table(newnames,nodepos,'VariableNames',{'Name','pos'});
         NT.graph = addnode(NT.graph,nodetab);
         newind = findnode(NT.graph,newnames);
         NT.graph.Nodes.pos(newind,:) = nodepos;         
                  
         if (~isempty(newedges))
             NT.graph = addedge(NT.graph,newedges(:,1),newedges(:,2));
             % new edge index
             newedgeind = findedge(NT.graph,newedges(:,1),newedges(:,2))';
             
             if (exist('edgepaths','var'))
                 setedgepaths = ~isempty(edgepaths);
             else
                 setedgepaths = false;
             end
             
             if (setedgepaths)
                 for ec = 1:length(newedgeind)
                     NT.graph.Edges.edgepath{newedgeind(ec)} = edgepaths{ec};
                 end
             else
                 NT.interpolateEdgePaths(2, newedgeind);
             end
         end
     end
     
     function interpolateEdgePaths(NT,npt,whichedges)
        % for a network graph, get edge paths as a linear interpolation
        % between the end points of each edge
        % npt = number of points to interpolate along each edge.
        % IMPORTANT: this resets edge lengths to be the distances between end points
        % (should fix at some point to have similar step sizes in
        % interpolation (different npt for different length edges))
                
        if (~exist('whichedges','var'))
            whichedges = 1:NT.nedge;
        end
        
        Nodes = NT.graph.Nodes;
        
        ipts = linspace(0,1,npt)'; % interpolation points
        for ec = whichedges
            n1 = NT.graph.Edges.EndNodes(ec,1); n2 = NT.graph.Edges.EndNodes(ec,2);
            p1 = Nodes.pos(n1,:); p2 = Nodes.pos(n2,:);                    
            
            NT.graph.Edges.edgepath{ec} = p1 + ipts*(p2-p1);            
        end
        
        % set up cumulative lengths
        NT.setCumEdgeLen(whichedges,true);
        
     end
     
     function breakEdge(NT,ectarget,breakfrac)
         % break up an edge (index ectarget) in the network, creating a new degree-2 node along it
         % new node is located at fraction breakfrac along the edge
        
         %ectarget = findedge(NT.graph,n1,n2);
         Edges = NT.graph.Edges; 
         prevnodes = Edges.EndNodes(ectarget,:);
         
         
        % get new node position
        shiftlen = Edges.edgelens(ectarget)*breakfrac;
        newpos = interp1(Edges.cumedgelen{ectarget},Edges.edgepath{ectarget},shiftlen);
        % last path index before this point
        ind = find(Edges.cumedgelen{ectarget}<shiftlen,1,'last');
        
        path1 = [Edges.edgepath{ectarget}(1:ind,:); newpos];
        path2 = [newpos; Edges.edgepath{ectarget}(ind+1:end,:)];
        
        % get rid of old edge
        NT.graph = rmedge(NT.graph,prevnodes(1),prevnodes(2));
        
        % add in new node
        newedgeind = NT.addNodes(newpos,{prevnodes},{path1,path2});   
        % recalculate lengths
        NT.setCumEdgeLen(newedgeind,true);
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
         
         Edges = NT.graph.Edges;
         nedge = size(Edges,1);
         
         if (~ismember('edgepath',Edges.Properties.VariableNames))
             error('need to set edge paths first')
         end
         
         distance = inf*ones(size(pts,1),1);
         for ec = 1:nedge
             curvexy = Edges.edgepath{ec};
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
    
     function redistributeEdgePaths(NT,dxwant)
         % reinterpolate edge paths to have node separations roughly dxwant
         
         Edges = NT.graph.Edges;
         
         for ec = 1:NT.graph.numedges
             npt = ceil(Edges.edgelens(ec)/dxwant) + 1;
             NT.reinterpolateEdgePath(ec,max(npt,2));
         end
     end
     
     function fixEdgePathEndPoints(NT,ec)
         % fix endpoints of an edge path to match up to the end-nodes of
         % this edge
         Edges = NT.graph.Edges;
         
         path = Edges.edgepath{ec};
         
         nc = findnode(NT.graph,Edges.EndNodes(ec,:));
         path(1,:) = NT.graph.Nodes.pos(nc(1),:);
         path(end,:) = NT.graph.Nodes.pos(nc(2),:);
         NT.graph.Edges.edgepath{ec} = path;
         
         NT.setCumEdgeLen(ec,true);
         
     end
     
     function reinterpolateEdgePath(NT,edgeind,npt)
         % redistribute points along an edge path          
         % also adjusts endpoints and cumlengths
         
         Edges = NT.graph.Edges;
         for ec = edgeind                          
             NT.graph.Edges.edgepath{ec} = splinePath(Edges.edgepath{ec},Edges.cumedgelen{ec},npt);                         
         end
         NT.setCumEdgeLen(edgeind,true);
         
     end
   
     function newec = mergeEdgesAtNode(NT,nodename)
         % merge edges across a degree 2 node, to make a single edge        
         
         G = NT.graph; Edges = G.Edges;
         
         nc = findnode(G,nodename);
         
         outE = outedges(G,nc);
         inE = inedges(G,nc);
         
         deg = length(outE)+length(inE);
         if (deg ~=2)
             error(sprintf('Can only merge at degree 2 nodes. This one has degree: %d', deg));
         end                  
             
         if (isempty(outE))
             % both edges in, new edge direction is arbitrary
             neighb = Edges.EndNodes(inE,1);                     
             newpath = [Edges.edgepath{inE(1)}; flipud(Edges.edgepath{inE(2)}(1:end-1,:))];             
         elseif (isempty(inE))
             % both edges out, new edge direction is arbitrary 
             neighb = Edges.EndNodes(outE,2);
             newpath = [flipud(Edges.edgepath{outE(1)}(2:end,:)); Edges.edgepath{outE(2)}];
         else 
             % one edge in, one out; keep direction
             neighb(1) = G.Edges.EndNodes(inE,1);
             neighb(2) = G.Edges.EndNodes(outE,2);
             newpath = [Edges.edgepath{inE}; Edges.edgepath{outE}(2:end,:)];
         end                  
        
         % add new edge
         newedgeind = NT.addEdges(neighb(1),neighb(2),{newpath});
         
         % remove old node and its edges
         keepnodes = true(NT.graph.numnodes,1);
         keepnodes(nc) = false;
         NT.graph = subgraph(NT.graph,keepnodes);
     end
     
     function [nodeIDs, nodepos] = nearestNode(NT,pos)
         % find the closest nodes (in euclidean distance) to the specified
         % positions
         % also returns the node positions
                  
         nodeIDs = dsearchn(NT.graph.Nodes.pos,pos);
         nodepos = NT.graph.Nodes.pos(nodeIDs,:);
     end
     
     function NTold = convert2OldNetwork(NT)
         % convert to an old-style network object
         
         % ordered list of node names
         nodenames = NT.graph.Nodes.Name;
         
         NTold = NetworkObj();
         NTold.nodepos = [NT.graph.Nodes.pos];        
         for pc = 1:2
             NTold.edgenodes(:,pc) = findnode(NT.graph,NT.graph.Edges.EndNodes(:,pc));
         end
         
         
         NTold.setupNetwork()
         for ec= 1:NTold.nedge
             NTold.edgepath{ec} = NT.graph.Edges.edgepath{ec};
             NTold.cumedgelen{ec} = NT.graph.Edges.cumedgelen{ec};
             NTold.edgelens(ec) = NT.graph.Edges.edgelens(ec);
         end
         
     end
     
     function getFromOldNetwork(NTg,NTold)
         % convert an old-style network object into a new one         
            NTg.setupNetwork(NTold.edgenodes,NTold.nodepos,NTold.edgepath');
     end         
         
     
    end
    
end