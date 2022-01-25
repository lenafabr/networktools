function sp = splinePath(path,cumlen,npt)
% given a path and cumulative arc lengths, get points along the
% interpolating spline
% path = nx2 coordinates of path points
% cumlen = cumulative arc lengths along path
% npt = number of points to interpolate to

    samppt = linspace(0,cumlen(end),npt);
             

    snakeIx = spline(cumlen,path(:,1),samppt);
    snakeIy = spline(cumlen,path(:,2),samppt);
    
    sp = [snakeIx; snakeIy]';
end