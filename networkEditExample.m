%% This example walks through the steps to go from a black-and-white
% segmented image (OR already skeletonized image) to a network structure, that can then be edited manually
% with a GUI
% followed by measurements of edge widths along the network structure
% To use, run through this matlab script, one section at a time (eg: by
% pressing Ctrl-Enter on each section)

%%
% set up your file paths (may require tweaking on a Windows/Mac system)

% make sure matlab knows about the gui directory
addpath('./gui')

% original image file, not segmented, will be used for visualization only
origimgfile = './examples/example2_dendriticNetwork.tif';
% single-frame image file for a black and white image (can be skeletonized
% already, or merely segmented)
bwimgfile = './examples/example2_dendriticNetworkSkeleton.tif';

% load the images
img = imread(origimgfile);
bwimg = imread(bwimgfile);

%% Extract a network object structure from the black and while image
% this step can be  slow for a large image and has not been fully optimized
NT= getNetworkFromBWImage(bwimg);
% make a backup copy of the network before you edit it
NT0 = copy(NT);
%% visualize the network object to make sure it looks ok

% plotting options (avoids drawing black lines on black background)
plotopt = struct('nodesize',20,'nodecolor',[1 0 0],'datatipindex',true)
plotopt.edgeplotopt = {'LineWidth',2,'Color','g'};

% superimpose network on image
% replace img with bwimg if you want to see how the network compares to the
% segmented image it was calculated from
imshow(img,[])
hold all
NT.plotNetwork(plotopt)
hold off

%% Use a graphical user interface (GUI) to edit the network (eg: remove or add nodes / edges)
% when you are done editing (or as you go), the NT object will contain the
% new network structure
networkEdit('NT',NT,'img',img,'plotopt',plotopt)


%% Now we will use a separate GUI to measure the edge widths and store the measurements
% You can measure an arbitrary number of separate widths on each edge
% [these two GUIs will be integrated into one eventually....]

% set up edgewidths array where widths will be stored
%% set dimensions of edgewidths

% run the gui for width measurement
% pink dots are locations of existing width measurements that are already
% saved in the object
[NT,saveEdgeVals] = setNetWidths(NT,origimgfile);
