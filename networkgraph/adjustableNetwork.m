classdef adjustableNetwork < handle
    % class defining interactive ROIs for adjusting a network
    
    properties
        lineColor='g' % Spline color property
        markerColor='r' % Spline control node marker color
        marker='s' % Spline control node marker shape
        edgemarker = 's'; % marker for edge control nodes
        edgemarkerColor = 'm'; 
        nptspline = 51; % number of points to use when drawing splines.
    end
    
     properties(SetAccess=private)
        parent % parent axes
        NT % associated network graph object        
        nodeH % handles to network nodes
        splineH % handles to edge splines
        edgeH % control points along edges
        
        nodeedges % keep track of edges connected to each node
     end
    
     methods
         function self = adjustableNetwork(NT,ax)
             % create a spline network object from a network object
             % (uses new NetworkGraphObj)
             
             self.NT = NT;
             %self.splineROIs = splineroi.empty();
             
             self.setNodeEdges()
             
             % set up axes
             if (exist('ax','var'))
                 self.parent = ax;
             else
                 self.parent = gca;
             end
             
             set(self.parent,'nextplot','add');
                 
             % plot edge splines
             self.splineH = gobjects(NT.graph.numedges,1);
             for ec = 1:self.NT.graph.numedges
                 self.plotSpline(ec);
             end
             
             % plot edge control points             
             self.edgeH = {};             
             for ec = 1:self.NT.graph.numedges                                  
                 path = NT.graph.Edges.edgepath{ec};
                 npt = size(path,1)-2;
                 self.edgeH{ec} = gobjects(npt,1);
                 for pc = 1:npt                    
                     self.plotEdgePoint(ec,pc);
                 end
             end
             
             % put point ROIs on network nodes     
             self.nodeH = gobjects(NT.graph.numnodes,1);
             for nc = 1:self.NT.graph.numnodes                 
                 self.plotNode(nc);
                %drawpoint('Position',Nodes.pos(nc,:),'Color','b');
             end

         end

         function delete(self)
             % clear all graphics handles
             delete(self.nodeH)
             for ec = 1:length(self.edgeH)
                delete(self.edgeH{ec})
             end
             delete(self.splineH)
         end
         
         function setNodeEdges(self)
             % keep track of edges connected to each node, to avoid
             % recalculating
             G = self.NT.graph;
             for nc = 1:G.numnodes
                 self.nodeedges{nc} = [outedges(G,nc); inedges(G,nc)];
             end
         end
         
         function self=plotNode(self,iNode)
            % FUNCTION plotNode
            % Plots one node and sets properties to default values
            if ~ishandle(self.nodeH(iNode))
                hNew = plot(self.NT.graph.Nodes.pos(iNode,1),self.NT.graph.Nodes.pos(iNode,2),'rs');
                self.nodeH(iNode) = hNew;
                set(self.nodeH(iNode),'markerfacecolor',self.markerColor,'markerEdgeColor',self.markerColor,'Marker',self.marker);
                set(self.nodeH(iNode),'buttondownfcn',{@adjustableNetwork.animator,self,'start','node',iNode});
            end
         end
         
          function self=plotEdgePoint(self,ec,iNode)
            % FUNCTION plotNode
            % Plots one node along an edge and sets properties
            % ec = which edge, iNode = which index along it
            if ~ishandle(self.edgeH{ec}(iNode))
                pt = self.NT.graph.Edges.edgepath{ec}(iNode+1,:);
                hNew = plot(pt(1),pt(2),'rs');
                self.edgeH{ec}(iNode) = hNew;
                set(self.edgeH{ec}(iNode),'markerfacecolor',self.edgemarkerColor,'markerEdgeColor',self.edgemarkerColor,'Marker',self.edgemarker);
                set(self.edgeH{ec}(iNode),'buttondownfcn',{@adjustableNetwork.animator,self,'start','edge',[ec iNode]});
            end
         end
        
         function self=plotSpline(self,ec)
            % FUNCTION plotSPline
            % Plots the spline          
            Edges= self.NT.graph.Edges;
            sPoints = splinePath(Edges.edgepath{ec},Edges.cumedgelen{ec},self.nptspline);
            %sPoints = [sPoints sPoints(:,1)];
            if ishandle(self.splineH(ec))
                set(self.splineH(ec),'xdata',sPoints(:,1),'ydata',sPoints(:,2))
            else
                self.splineH(ec)=plot(sPoints(:,1),sPoints(:,2),'color',self.lineColor,'lineStyle','-','LineWidth',2);
                set(self.splineH(ec),'hittest','off');
            end            
        end
     end
     
     methods(Static)
        function animator(src,eventdata,self,action,type,idx) %#ok<INUSL>
            % FUNCTION animator
            % Adds drag and drop functionality to spline nodes
            % idx = edge index and node index within the edge            
            switch(action)
                case 'start'
                    set(gcbf,'WindowButtonMotionFcn',{@adjustableNetwork.animator,self,'move',type,idx});
                    set(gcbf,'WindowButtonUpFcn',{@adjustableNetwork.animator,self,'stop',type,idx});
                    %self.doSplineChangedNotification = false;
                case 'move'
                    currPt=get(gca,'CurrentPoint');
                    set(gco,'XData',currPt(1,1));
                    set(gco,'YData',currPt(1,2));
                    if (strcmp(type,'node'))
                        self.NT.graph.Nodes.pos(idx,:) = currPt(1,1:2);
                        
                        edges = self.nodeedges{idx}';
                        for ec = edges
                            % update and replot edge
                            %npt = size(self.NT.graph.Edges.edgepath{ec},1); % number of edge control points
                            self.NT.fixEdgePathEndPoints(ec)
                            self.plotSpline(ec);
                        end
                    elseif (strcmp(type,'edge'))                        
                        self.NT.graph.Edges.edgepath{idx(1)}(idx(2)+1,:) = currPt(1,1:2);                        
                        % reset cumedgelengths
                        self.NT.setCumEdgeLen(idx(1));
                        self.plotSpline(idx(1));    
                    else
                        error('type must be node or edge')
                    end
                   % self.splineNodes(self.hPoint==gco,:)=currPt(1,1:2);
                   % self.plotSpline;                 
                case 'stop'
                    set(gcbf,'WindowButtonMotionFcn','');
                    set(gcbf,'windowButtonUpFcn','');
                    %self.doSplineChangedNotification = true;
                    %notify(self,'splineChanged');
            end%switch
        end
    end
      
end