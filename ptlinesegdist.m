function [dist,t,ptnear] = ptlinesegdist(pt,linept1,linept2)
%% distance from point to line segment
% row vectors assumed
% NOT optimized for speed

v = linept2-linept1;
nv = norm(v);

d = pt-linept1;
nd = norm(d);
rho = d*v'/nd/nv;

% fraction along segment
t = nd*rho/nv;
if (t<0)
    % nearest distance is to first endpoint
    dist = nd;
    ptnear = linept1;
elseif (t>1)
    % nearest distance is to last endpoint
    dist = norm(pt-linept2);
    ptnear = linept2;
else
    % nearest distance to line
    dist = nd*sqrt(1-rho^2);
    ptnear = t*v+linept1;
end