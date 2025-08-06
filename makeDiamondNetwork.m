%% Create a diamond cubic lattice spherical shell network
% This script builds a 3D diamond lattice, truncates it to a spherical shell,
% adds radial tethers (inward and outward), wraps into NetworkObj, and writes a .net file.

% Parameters
a = 2;                          % Diamond lattice constant
sphereRadius = 8;              % Outer radius of the spherical shell
coreRadius = 4;                % Inner hollow core radius
outerTetherLength = 1;         % Length of outward tethers
innerTetherLength = 0.5;       % Length of inward tethers

% Diamond basis (FCC + tetrahedral shift)
fcc_basis = [0, 0, 0;
             0, 0.5, 0.5;
             0.5, 0, 0.5;
             0.5, 0.5, 0];

shift = [0.25, 0.25, 0.25];
diamond_basis = [fcc_basis; fcc_basis + shift];  % 8 atoms per cell
diamond_basis = unique(diamond_basis, 'rows') * a;

% Lattice vectors
a1 = a * [1, 0, 0];
a2 = a * [0, 1, 0];
a3 = a * [0, 0, 1];

% Generate nodes inside spherical shell
ncell = ceil((sphereRadius + 2) / a);
nodes = [];
for ix = -ncell:ncell
    for iy = -ncell:ncell
        for iz = -ncell:ncell
            offset = ix*a1 + iy*a2 + iz*a3;
            for b = 1:size(diamond_basis,1)
                pos = offset + diamond_basis(b,:);
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

% Build edges: nearest neighbors in diamond are at a*sqrt(3)/4
target_dist = sqrt(3)/4 * a;
cutoff = 1.01 * target_dist;
edges = [];
for i = 1:nNodes
    for j = i+1:nNodes
        dist = norm(nodes(i,:) - nodes(j,:));
        if abs(dist - target_dist) < 1e-3
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
title('Diamond Shell with Radial Tethers (FN = red, R1 = green)');
view(3); grid on;

%% Wrap into NetworkObj
NT = NetworkObj();
NT.nodepos = nodes;
NT.edgenodes = edges;
NT.nodelabels = labels;
NT.setupNetwork();

% Optional: Visual check
figure;
NT.plotNetwork(struct('labels', 0));
title('NetworkObj NT: Diamond Lattice');

%% Store it in the .net file
output_filename = '/home/yuz261/Research/ER_Ca_refill_sims/diamond/diamond.net';
fid = fopen(output_filename, 'w');

% Write node section
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

% Write edge section
fprintf(fid, '\n# list of edges (index node1 node2 length)\n');
for i = 1:NT.nedge
    n1 = NT.edgenodes(i,1);
    n2 = NT.edgenodes(i,2);
    len = NT.edgelens(i);
    fprintf(fid, 'EDGE %d %d %d %0.10f\n', i, n1, n2, len);
end

fclose(fid);
fprintf('✅ Diamond lattice network saved to: %s\n', output_filename);

%% Check the degree of nodes for debug
% Ensure setup was done
NT.setupNetwork();  % <- make sure this is called before!

% Get node degrees
nodeDegrees = NT.degrees;

% Count frequency of each degree
maxDeg = max(nodeDegrees);
degBins = 0:maxDeg;
degCounts = histcounts(nodeDegrees, [degBins, maxDeg+1]);  % bin edges must be 1 longer than bin count

% Plot
figure;
bar(degBins, degCounts, 'FaceColor', [0.2 0.5 0.8]);
xlabel('Node Degree');
ylabel('Number of Nodes');
title('Degree Distribution of Nodes in the Network');
grid on;

