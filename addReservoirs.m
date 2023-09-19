function [NTnew, mapold2new] = addReservoirs(NT,resvbound,minedgelen)

% add a reservoir to a network
% bound = nx2 array of points defining the reservoir polygon
% any nodes in the reservoir get removed
% new nodes formed at the half-way point of edges crossing the reservoir
% boundary
% Note that this will give a reservoir of a slightly different shape than
% input, *but* will avoid formation of any super-short edges

% bound can be a cell array with multiple regions

% minedgelen: when breaking up edges at the reservoir boundary,
% do not make any shorter than this (or half of previous edge length if
% smaller)

if ~iscell(resvbound)
    resvbound = {resvbound};
end

ct = 0;
newnodes = []; newcon = {}; sheet2node = {}; newedgepath={};
dokeepnode = true(1,NT.nnode);

for sc = 1:length(resvbound)
    sheet2node{sc} = [];
    
    bound = resvbound{sc};
    boundw = [bound; bound(1,:)]; % wrap around
    
    % find nodes within the reservoir boundary
    nodeisin = inpolygon(NT.nodepos(:,1),NT.nodepos(:,2), bound(:,1),bound(:,2));        
    nodelist = find(nodeisin);
    
    for ic = 1:length(nodelist)
        nc = nodelist(ic);
        
        dokeepnode(nc) = false; % will drop nodes inside the sheet
        
        p1 = NT.nodepos(nc,:);
        neighb = NT.nodenodes(nc,1:NT.degrees(nc));

        for ic2 = 1:length(neighb)
            n2 = neighb(ic2);
            % only look at neighbors outside the reservoir
            if (nodeisin(n2)); continue; end
            p2 = NT.nodepos(n2,:);

            % edge going to this neighbor
            ec = NT.nodeedges(nc,ic2);
            minlenuse = min(minedgelen, NT.edgelens(ec)/2);
            
            path = NT.edgepath{ec};
            
            % intersection of edge path with the boundary          
            [xint,yint,ii] = polyxpoly(NT.edgepath{ec}(:,1),NT.edgepath{ec}(:,2),boundw(:,1),boundw(:,2));
            xint = xint(1); yint = yint(1); ii = ii(1);
            pint = [xint yint];
            
            ct = ct+1; % count new nodes formed            
            
            breaklen = NT.cumedgelen{ec}(ii) + norm(pint-path(ii,:));                                               
            if (NT.edgenodes(ec,1)==n2) % inward going edge
                if (breaklen<minlenuse) % break at minimum length along edge path
                    shiftlen = minlenuse;
                    newpos = interp1(NT.cumedgelen{ec},[path (1:size(path,1))'],shiftlen);
                    pint = newpos(:,1:2); indprev = floor(newpos(3));                    
                else
                    pint = [xint yint]; indprev = ii;                    
                end
                
                % distance from previous pt
                fromprev = norm(pint - path(indprev,:));                
                % path from intersection to outer point
                if (indprev>1 & fromprev<1e-6)
                    % breaking right above control point
                    newedgepath{ct}= [pint; flipud(path(1:indprev-1,:))];
                else
                    newedgepath{ct}= [pint; flipud(path(1:indprev,:))];
                end
            else % outward going edge                
                 if (breaklen > NT.cumedgelen{ec}(end) - minlenuse) % break half-way along edge path
                    shiftlen = NT.cumedgelen{ec}(end) - minlenuse;
                    newpos = interp1(NT.cumedgelen{ec},[path (1:size(path,1))'],shiftlen);
                    pint = newpos(:,1:2); indprev = floor(newpos(3));
                 else
                     pint = [xint yint]; indprev = ii;
                 end
                 
                 % distance from next pt
                 fromnext = norm(pint - path(indprev+1,:));
                 
                 % path from intersection to outer point
                if (indprev<size(path,1)-1 & fromnext<1e-6)
                    % breaking right below control point
                    newedgepath{ct}= [pint; path(indprev+2:end,:)];
                else
                    newedgepath{ct}= [pint; path(indprev+1:end,:)];
                end
                
            end
                    
            
            % create a new node and new connection to it            
            newnodes(ct,:) = pint;           
            newcon{ct} = [n2];       
           
            sheet2node{sc} = [sheet2node{sc} ct+NT.nnode];                        
        end
    end
end

%% add new nodes to the network
nnprev = NT.nnode;
neprev = NT.nedge;

NTnew = copy(NT);
NTnew.addNodes(newnodes,newcon);
nn = NTnew.nnode;
ne = NTnew.nedge;

for cc = 1:length(newedgepath)   
    NTnew.edgepath{neprev+cc} = newedgepath{cc};
end
NTnew.setCumEdgeLen(1:NTnew.nedge,true);

allkeepnode = true(1,NTnew.nnode);
allkeepnode(find(~dokeepnode)) = false;

%% remove nodes within the sheets
mapold2new = NTnew.keepNodes(find(allkeepnode));
%
for sc = 1:length(sheet2node)
    sheet2node{sc} = mapold2new(sheet2node{sc});
    
    % label nodes belonging to reservoir
    for ic = 1:length(sheet2node{sc})
        nc = sheet2node{sc}(ic);
        NTnew.nodelabels{nc} = sprintf('R%d',sc);
    end
end


