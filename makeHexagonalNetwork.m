%% Create a hexagonal close-packed (HCP) spherical shell network
% This script builds a 3D HCP lattice, truncates it to form a spherical shell,
% and adds radial tethers pointing outward and inward.

% Parameters
a = 2;                          % Lattice spacing between nearest neighbors
c = sqrt(8/3)*a;                % Vertical spacing for ideal HCP lattice (c/a ratio ≈ 1.633)

sphereRadius = 8;              % Outer radius of the spherical shell (defines the whole structure size)
coreRadius = 4;                % Inner hollow radius (nodes inside this sphere are removed)

outerTetherLength = 1;         % Length of tethers sticking radially outward from outer surface nodes
innerTetherLength = 0.5;       % Length of tethers pointing inward from inner surface nodes toward the center


% HCP basis & vectors
basis = [0, 0, 0;
         0.5*a, sqrt(3)/6*a, 0.5*c];
a1 = [a, 0, 0]; a2 = [0.5*a, sqrt(3)/2*a, 0]; a3 = [0, 0, c];

% Generate nodes inside spherical shell
ncell = ceil((sphereRadius + 2) / a);
nodes = [];
for ix = -ncell:ncell
    for iy = -ncell:ncell
        for iz = -ncell:ncell
            offset = ix*a1 + iy*a2 + iz*a3;
            for b = 1:size(basis,1)
                pos = offset + basis(b,:);
                r = norm(pos);
                if coreRadius < r && r <= sphereRadius
                    nodes(end+1,:) = pos;
                end
            end
        end
    end
end
nodes = unique(round(nodes, 6), 'rows');
nNodes = size(nodes,1);
labels = repmat({''}, size(nodes,1), 1);

% Build normal edges
edges = [];
cutoff = 1.01 * a;
for i = 1:nNodes
    for j = i+1:nNodes
        dist = norm(nodes(i,:) - nodes(j,:));
        if abs(dist - a) < 1e-3 || abs(dist - c) < 1e-3
            edges(end+1,:) = [i, j];
        end
    end
end

% Surface node detection
radii = vecnorm(nodes, 2, 2);
surfaceOuter = find(abs(radii - sphereRadius) < a/2);
surfaceInner = find(abs(radii - coreRadius) < a/2);

% Add outward tethers ("FN")
for i = 1:length(surfaceOuter)
    idx = surfaceOuter(i);
    pos = nodes(idx, :);
    dir = pos / norm(pos);
    newNode = pos + outerTetherLength * dir;
    nodes(end+1,:) = newNode;
    edges(end+1,:) = [idx, size(nodes,1)];
    labels{end+1} = 'FN';
end

% Add inward tethers ("R1")
for i = 1:length(surfaceInner)
    idx = surfaceInner(i);
    pos = nodes(idx, :);
    dir = -pos / norm(pos);
    newNode = pos + innerTetherLength * dir;
    nodes(end+1,:) = newNode;
    edges(end+1,:) = [idx, size(nodes,1)];
    labels{end+1} = 'R1';
end

% Assign colors based on labels
colors = zeros(size(nodes,1), 3); % default blue
for i = 1:length(labels)
    if strcmp(labels{i}, 'FN')
        colors(i,:) = [1, 0, 0]; % red
    elseif strcmp(labels{i}, 'R1')
        colors(i,:) = [0, 0.6, 0]; % green
    else
        colors(i,:) = [0, 0, 1]; % blue
    end
end

% Plot
figure; hold on;
for e = 1:size(edges,1)
    p1 = nodes(edges(e,1), :);
    p2 = nodes(edges(e,2), :);
    plot3([p1(1), p2(1)], [p1(2), p2(2)], [p1(3), p2(3)], 'k-');
end
scatter3(nodes(:,1), nodes(:,2), nodes(:,3), 20, colors, 'filled');
axis equal
xlabel('X'); ylabel('Y'); zlabel('Z');
title('HCP Shell with Radial Tethers (FN = red, R1 = green)');
view(3); grid on;


%% Wrap into NetworkObj
NT = NetworkObj();
NT.nodepos = nodes;        % node positions
NT.edgenodes = edges;      % edge list
NT.nodelabels = labels;    % <<< assign labels before setup
NT.setupNetwork();         % computes degrees, edge lengths, etc.

% View with labeled node colors
figure;
NT.plotNetwork(struct('labels', 0));  % show node indices if needed
title('NetworkObj NT with Node Labels');


%% Store it in the .net file
output_filename = '/home/yuz261/Research/ER_Ca_refill_sims/hexa/hexa.net';
fid = fopen(output_filename, 'w');

% --- Write node section ---
fprintf(fid, '# list of node indices and xyz positions and values\n');
for i = 1:NT.nnode
    pos = NT.nodepos(i,:);
    if length(pos) < 3
        pos(3) = 0;
    end
    label = '';
    if isprop(NT, 'nodelabels') && ~isempty(NT.nodelabels) && ischar(NT.nodelabels{i}) && ~strcmp(NT.nodelabels{i}, '')
        label = NT.nodelabels{i};
        fprintf(fid, 'NODE %d %0.10f %0.10f %0.10f 0 %s\n', i, pos(1), pos(2), pos(3), label);
    else
        fprintf(fid, 'NODE %d %0.10f %0.10f %0.10f 0\n', i, pos(1), pos(2), pos(3));
    end
end

% --- Write edge section ---
fprintf(fid, '\n# list of edges (index node1 node2 length)\n');
for i = 1:NT.nedge
    n1 = NT.edgenodes(i,1);
    n2 = NT.edgenodes(i,2);
    len = NT.edgelens(i);
    fprintf(fid, 'EDGE %d %d %d %0.10f\n', i, n1, n2, len);
end

fclose(fid);
fprintf('✅ Network written to %s\n', output_filename);


