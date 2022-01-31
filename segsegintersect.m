function [t1,t2,intpt,doint] = segsegintersect(seg1, seg2)
% find the intersection of 2 intersecting line segments
% seg(:,i) are the coords of endpoint i of the segment
% return t1 = segment 1 parameter at intersection (btwn 0 and 1)
% t2 = segment 2 parameter at intersection
% intpt the coordinates of intersection
% if t1 or t2 outside of 0,1 then do not intersect
% doint is 1 if they do intersect and 0 otherwise


intpt = [0;0];

ds1 = seg1(:,2)-seg1(:,1);
ds2 = seg2(:,2)-seg2(:,1);

mat = [ds1(1) -ds2(1); ds1(2) -ds2(2)];
vec = seg2(:,1) - seg1(:,1);

tmp = mat\vec;

t1 = tmp(1); t2 = tmp(2);

doint = (t1>=0 & t1<=1 & t2>=0 & t2<=1);
if (doint)
    intpt = seg2(:,1) + t2*ds2;
end

end