function [distmin, breakfracmin] = setNodeNearestPoint(NT,point,minedgelen)
% modify a network to create a node along the curved edges
% at the point closest to a given probe point

npt = 100;
distmin = inf;
for ec = 1:NT.nedge
    xp = linspace(0,NT.edgelens(ec),npt);
    ptint = interp1(NT.cumedgelen{ec},NT.edgepath{ec},xp);
    distint = sqrt(sum((ptint-point).^2,2));
    [tmp,b] = min(distint);
    breakfrac = xp(b)/NT.edgelens(ec);
    
    if tmp<distmin
        distmin = tmp;
        ecmin = ec;
        breakfracmin = breakfrac;
    end
end

% only break the edge if this does not generate very short edges
elen = NT.edgelens(ecmin);
if (breakfracmin*elen<minedgelen | (1-breakfracmin)*elen<minedgelen)
    % do not break network
    display('cannot break here, edge is too short')
    return
else
    % create a new point along the edge
    NT.breakEdge(ecmin,breakfracmin,true);
end

end