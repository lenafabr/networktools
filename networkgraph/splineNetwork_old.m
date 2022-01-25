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
                 self.splineROIs(ec)= splineroi(self.parent);
                 self.splineROIs(ec).addNode(pts);
             end
         end
     end
     
     methods(Access=private)
        function self=plotNode(self,iNode)
            % FUNCTION plotNode
            % Plots one node and sets properties to default values

            if ~ishandle(self.hPoint(iNode))
                hNew = plot(self.splineNodes(iNode,1),self.splineNodes(iNode,2),'rs');
                self.hPoint(iNode) = hNew;
                set(self.hPoint(iNode),'markerfacecolor',self.markerColor,'markerEdgeColor',self.markerColor,'Marker',self.marker);
                set(self.hPoint(iNode),'buttondownfcn',{@splineroi.animator,self,'start'});
            end
        end
        
        function self=plotSpline(self)
            % FUNCTION plotSPline
            % Plots the spline
            %sPoints = interpolateSnake(self, 101);
            %sPoints = sPoints';
            sPoints = self.calcSpline;
            %sPoints = [sPoints sPoints(:,1)];
            if ishandle(self.hSpline)
                set(self.hSpline,'xdata',sPoints(1,:),'ydata',sPoints(2,:))
            else
                self.hSpline=plot(sPoints(1,:),sPoints(2,:),'color',self.lineColor,'lineStyle','-');
                set(self.hSpline,'hittest','off');
            end            
        end
    end
    
    methods(Static)
        function animator(src,eventdata,self,action) %#ok<INUSL>
            % FUNCTION animator
            % Adds drag and drop functionality to spline nodes
            switch(action)
                case 'start'
                    set(gcbf,'WindowButtonMotionFcn',{@splineroi.animator,self,'move'});
                    set(gcbf,'WindowButtonUpFcn',{@splineroi.animator,self,'stop'});
                    self.doSplineChangedNotification = false;
                case 'move'
                    currPt=get(gca,'CurrentPoint');
                    set(gco,'XData',currPt(1,1));
                    set(gco,'YData',currPt(1,2));
                    self.splineNodes(self.hPoint==gco,:)=currPt(1,1:2);
                    self.plotSpline;
                case 'stop'
                    set(gcbf,'WindowButtonMotionFcn','');
                    set(gcbf,'windowButtonUpFcn','');
                    self.doSplineChangedNotification = true;
                    notify(self,'splineChanged');
            end%switch
        end
    end
    events
        splineChanged
    end
end