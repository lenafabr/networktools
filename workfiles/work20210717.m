%% load segmented image
bwimgfile = '~/UCSD/data/Obara/Sec61-TA/Halo-Sec61-TA-100Hz_004b_2_Mask.tif';
bwimg = imread(bwimgfile,1)';
imshow(bwimg)

%% extract network structure

[NT,skelimage,opt] = getNetworkFromBWImage(bwimg,struct('dodisplay',2));

%% show network
NT.plotNetwork()

%% load original image
imgfile = '~/UCSD/data/Obara/Sec61-TA/Halo-Sec61-TA-100Hz_004b_2_maxS2N.tif';
img = imread(imgfile,1);

% renormalize
img = img/max(img(:));

%% superimpose network and image
imshow(img,[0,0.5])
hold all
plotopt = struct('nodesize',20,'nodecolor',[1 0 0])
plotopt.edgeplotopt = {'LineWidth',2,'Color','g'};
NT.plotNetwork(plotopt)
hold off

%% make example file
save('workfiles/example_network_Halo_Sec61_4b.mat','imgfilename','NT','img','bwimg','extractopt','plotopt')