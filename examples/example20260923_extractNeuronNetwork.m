%% access code path (adjust to whereever you have the networktools code)
addpath(genpath('../../networktools'))

%% load in images
% directory containing images
imagedirname = '../data/WildongerGroup/';

% load in original image
origimg = imread([imagedirname 'dendrites_solo.tif']);

% load in skeletonized binary image
bwimg0 = imread([imagedirname 'dendrites_mask.tif']);

% slight erosion of bw image
% this avoids fork-like artifacts at the branch ends
se = strel("disk",1); % 1-px disc erosion element
bwimg = imerode(bwimg0,se);
imshowpair(origimg,bwimg)

%% Extract network
% do not keep just the largest connected component. Allow disconnected
% networks
options = struct('keepconncomp',false);
[NT,skelimage,opt] = getNetworkFromBWImage(bwimg,options);

%% Display image and overlaid network
imshow(origimg)
plotopt = struct('nodesize',20,'nodecolor',[1 0 0],'datatipindex',true)
plotopt.edgeplotopt = {'LineWidth',2,'Color','g'};

hold all
NT.plotNetwork(plotopt)
hold off

%% Use a graphical user interface (GUI) to edit the network (eg: remove or add nodes / edges)
% when you are done editing (or as you go), hit Update Network the NT object will contain the
% new network structure
% NOTE: for use with mito distribution code, the network needs to have a
% well-specified trunk edge. So you need to make an extra node + trunk in the
% cell body to serve this purpose
% Also, use the gui to get rid of any degree 4 nodes
% When all done, set the parent root node
plotopt = struct('nodesize',20,'nodecolor',[1 0 0],'datatipindex',true)
plotopt.edgeplotopt = {'LineWidth',2,'Color','g'};
app = networkEdit_App('NT',NT,'img',origimg,'plotopt',plotopt);

%%
NT = NetworkObj('../data/WildongerGroup/WildongerGroupDAneuron_20260924.net',struct('dim',2));
imshow(origimg,[])
hold all
NT.plotNetwork(struct('datatipindex',1,'plotoverimage',1))
hold off

%% Clean up after the GUI
NT.rootnode = 185; % if reloading from .net file, will need to set root manually
NT.edgevals = []; % edgevals has no meaning here
NT.mergeAllEdgePaths(); % get rid of degree 2 nodes
NT.keepLargestConnComp(); % get rid of disconnected pieces

% make sure tree is properly directed from root node
directedTreeEdges(NT,NT.rootnode)
%% View the edited network. Star indicates root node
figure
imshow(origimg)

hold all
NT.plotNetwork(plotopt)
hold all
plot(NT.nodepos(NT.rootnode,1),NT.nodepos(NT.rootnode,2),'m*','MarkerSize',10,'LineWidth',2)
hold off


%% Save network for later (adjust to whatever directory you want to save in)
savedirname = '../data/WildongerGroup';
save([savedirname 'DAneuron_20261005.mat'],'NT')
NT.outputNetwork([savedirname 'DAneuron_20261005.net'],struct('WRITEPATHS',true));

%% Plot histogram of edge lengths, as an example statistic
figure
histogram(NT.edgelens,20)
xlabel('edge length (px)')
ylabel('count')