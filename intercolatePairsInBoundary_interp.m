function [regbound2, boundnodeind2] = intercolatePairsInBoundary_interp(regbound,boundnodeind,tubeconwidth,minedgelen)
% CLEAN THIS UP LATER TO USE PROPER SPLINE INTERPOLATION AROUND BOUNDARY
error('DO not use. Not written yet!')
% for a boundary polygon with some points marked as belonging to network
% nodes, make a boundary with control points surrounding those node
% connections
% will try to avoid boundary edges with less than minimal length
% but might happen anyway if tube connections are too close

regboundw =[regbound; regbound(1,:)]; % wrap around
boundnodeindw = [boundnodeind boundnodeind(1)];

% arc length parameterization
param0 = arclenparam(regboundw);

% pick out points at appropriate arc length on either side of node
% connections
nind = find(boundnodeind>0);
for pc = nind

end



%%

regbound2 = []; boundnodeind2 = [];
for pc = 1:size(newregbound,1)
    if (boundnodeind(pc)==0)
        % ordinary boundary point
        regbound2(end+1,:) = newregbound(pc,:);
        if (length(boundnodeind2)<size(regbound2,1))
            % only update if not already taken over by prior node
            boundnodeind2(end+1) = 0;
        end
    else
        % node connection point
        if (pc == 1)
            v1 = newregboundw(end,:) - newregboundw(pc,:);
        else
            v1 = newregboundw(pc-1,:) - newregboundw(pc,:);
        end
        d1 = norm(v1);
        if (d1>minedgelen) % insert new point
            regbound2(end+1,:) = newregboundw(pc,:) + v1/d1*tubeconwidth;
            boundnodeind2(end+1) = boundnodeind(pc);
        elseif (boundnodeind2(end)>0)
            % the other point is already surrounding a tube, have to insert anyway
            regbound2(end+1,:) = newregboundw(pc,:) + v1/3;
            boundnodeind2(end+1) = boundnodeind(pc);
        else
            boundnodeind2(end) = boundnodeind(pc);
        end

        v2 = newregboundw(pc+1,:) - newregboundw(pc,:);
        d2 = norm(v2);

        if (d2>minedgelen) % insert new point
            regbound2(end+1,:) = newregboundw(pc,:) + v2/d2*tubeconwidth;
            boundnodeind2(end+1) = boundnodeind(pc);
        elseif (boundnodeindw(pc+1)>0)
            % the other point is already surrounding a tube, have to insert anyway
            regbound2(end+1,:) = newregboundw(pc,:) + v2/3;
            boundnodeind2(end+1) = boundnodeind(pc);
        else
            if (pc == size(newregbound,1))
                boundnodeind2(1) = boundnodeind(pc);
            else
                boundnodeind2(end+1) = boundnodeind(pc);
            end
        end
    end
end
end