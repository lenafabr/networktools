classdef splineNetwork < handle
    % class defining multiple interactive splineroi objects linked into a
    % network
    
     properties(SetAccess=private)
        parent
        splineROIs % spline ROIs for each edge
        NT % associated network graph object        
        
     end
    
     methods
         function self = splineNetwork(NT,ax)
             % create a spline network object from a network object
             % (uses new NetworkGraphObj)
             
             self.NT = NT;
             self.splineROIs = splineroi.empty();
             
             % set up axes
             if (exist('ax','var'))
                 self.parent = ax;
             else
                 self.parent = gca;
             end
                 
             % put point ROIs on network nodes
%              Nodes = NT.graph.Nodes;
%              for nc = 1:NT.graph.numnodes
%                 self.netnodeH(nc) = drawpoint('Position',Nodes.pos(nc,:),'Color','b');
%              end
             
           %  make spline for each edgepath
             Edges = NT.graph.Edges;
             for ec = 1:NT.graph.numedges
                 pts = Edges.edgepath{ec};
                 edgespline = splineroi(self.parent,self,ec);
                 edgespline.addNode(pts);
                 self.splineROIs(ec) = edgespline;
             end
         end
     end
     
   
    events
        splineChanged
    end
end