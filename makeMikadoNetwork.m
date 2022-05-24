function NT = makeMikadoNetwork(Nrod,Lrod,Lspace)
% construct a mikado network by placing Nrod rods of length Lrod, randomly 
% oriented with their centers within a box of size Lspace. Intersections of
% rods determine nodes of mikado network. 

% In this version, do not keep the degree one nodes at the tips of edges,
% could add the functionality
%%

cpts=rand(Nrod,2)*Lspace; % center points of the rods
dirs=rand(Nrod,2)-.5;
dirs=Lrod.*dirs./sqrt(sum(dirs.^2,2)); % vectors of length Lrod determining direction of each rod

epts = [cpts - (dirs/2) cpts+dirs/2];

%% Make Mikado network
NT=NetworkObj();

A=zeros(size(epts,1));
B=zeros(size(epts,1));

ct=0; % running count of intersections
for ec=1:size(epts,1)-1
    for ec2=(ec+1):size(epts,1)
%         if (ec==ec2)
%             continue
%         end


        [xi,yi] = polyxpoly(epts(ec,[1 3]), epts(ec,[2 4]),epts(ec2,[1 3]) ,epts(ec2,[2 4]));
        if (~isempty(xi))
            A(ec,ec2) = sqrt(sum((epts(ec,[1 2]) - [xi yi]).^2,2))/Lrod; % fraction along edge to place node
            A(ec2,ec) = sqrt(sum((epts(ec2,[1 2]) - [xi yi]).^2,2))/Lrod; % fraction along edge to place node

            ct=ct+1;
            NT.nodepos(end+1,:) = [xi,yi];
            B(ec,ec2) = ct; % node number
        end
    end
end
B=B+triu(B,1)';

%% set edgenodes using connection matrix from above
for ec=1:(size(epts,1))

    if (nnz(A(ec,:))<=1)
        % only d1 or no nodes
        continue
    else
        % find 
        inds1 = find(A(ec,:)>0);
        [~,Is] = sort(A(ec,inds1));
        
        for ii=1:(length(Is)-1)
            n1=B(ec,inds1(Is(ii)));
            n2=B(ec,inds1(Is(ii+1)));

            NT.edgenodes(end+1,:) = [n1 n2];
        end

    end

end

%% setup network
NT.setupNetwork()
NT.setEdgeLens;

end
