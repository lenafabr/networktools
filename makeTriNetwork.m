function NT = makeTriNetwork(N,celldiam)
%% make a triangular network (degree 6 nodes)
% with unit edge lengths
% makes a hexagonal shaped network with center->vertex radius =  N edges

elen = celldiam/2/N;

% unit vectors
v1 = [1 0]*elen;
v2 = [0.5,sqrt(3)/2]*elen;


nlist = -N:N

ct=0;
nodepos = zeros(N^2,2);
ncoords = zeros(N^2,2);
for n1 = nlist
    if (n1>0)
        n2list = -N:(N-n1);
    else
        n2list = -(N+n1):N;
    end    
    for n2 = n2list
        ct = ct+1;
        % position of each node and coordinates as unit vector multiples
        nodepos(ct,:) = n1*v1+n2*v2;
        ncoords(ct,:) = [n1,n2];        
    end
end

nnode = ct;
%plot(nodepos(:,1),nodepos(:,2),'.-')
%axis equal

% connect up neighbor nodes
ect = 0;
edgenodes = [];
for n1 = 1:nnode
    for n2 = n1+1:nnode        
        coorddiff = ncoords(n1,:)-ncoords(n2,:);
        if (max(abs(coorddiff)) == 1 & (any(coorddiff==0) | sum(coorddiff)==0 ))
            ect = ect+1;
            edgenodes(ect,:) = [n1,n2];
        end
    end
end

% make network object
NT = NetworkObj();
NT.nodepos = nodepos; 
NT.edgenodes = edgenodes;
NT.setupNetwork();
NT.interpolateEdgePaths(2);
NT.setCumEdgeLen();
NT.plotNetwork()
end