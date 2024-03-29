function widthinfo = makeEdgeWidth(NT,ec,w,d)
    % width info line defining the width of a given edge
    % default if d is not specified: assume measurement is made at center of the edge    
    % cumulative lengths array must be previously set
    % does *not* alter NT object directly
    % widthinfo is a 1x6 vector listing: width, position along edge where
    % the width was calculated, x and y for first endpoint of width segment
    % then x and y for second endpoint of width segment
    
    path = NT.edgepath{ec};
    cum = NT.cumedgelen{ec};   

    % position along edge
    if (~exist('d','var'))
        d = NT.edgelens(ec)/2;
    end
    edgepos = interp1(cum,0:(length(cum)-1),d);
    ind = floor(edgepos)+1;
    frac = mod(edgepos,1);
    xypos = path(ind,:)*(1-frac) + path(ind+1,:)*frac;    

    % get perpendicular unit vector
    dvec = path(ind+1,:) - path(ind,:);
    dvec = dvec/norm(dvec);
    perpvec = [-dvec(2),dvec(1)];

    xy1 = xypos + perpvec*w/2;
    xy2 = xypos - perpvec*w/2;

    widthinfo = [w,d,xy1 xy2];    
end