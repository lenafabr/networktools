% play around with dendritic tree

%% load data
origimgfile =  '../dendriticmito/imgData/MAX_6.211020-HS_brain9_25min.tif';
img = imread(origimgfile);
img = img(:,:,1);

load('../dendriticmito/imgData/6.211020-HS_brain9_25min_network.mat','NT')

%%
load('~/proj/dendriticmito/results/workspace_filteredNT_20210817.mat')
imgfilename = '~/proj/dendriticmito/data/MCFO-HSN-A61-2.tif';
img = imread(imgfilename);



%% visualize

plotopt = struct('nodesize',20,'nodecolor',[1 0 0],'datatipindex',true)
plotopt.edgeplotopt = {'LineWidth',1,'Color','g'};
%%
% superimpose network on image
% replace img with bwimg if you want to see how the network compares to the
% segmented image it was calculated from
%imshow(bwimg,[])
%hold all
figure(2)
imshow(img,[])
%imshow(img(:,:,:))
hold all
NT.plotNetwork(plotopt)
hold off

%% save workspace
save('tmp.mat')
%% GUI to run network edit code - save existing NT workspace first!
% Need to load mat file before editing
networkEdit('NT',NT,'img',img,'plotopt',plotopt)

% -------- New Network Objects --------------

%% convert to new network
NTg = NetworkGraphObj();
NTg.getFromOldNetwork(NT);
NTg0 = copy(NTg);

%% visualize new network
plotopt = struct('nodesize',5,'nodecolor',[1 0 0],'datatipindex',true)
plotopt.edgeplotopt = {'LineWidth',1,'Color','g'};

imshow(img,[])
%imshow(img(:,:,:))
hold all
NTg.plotNetwork(plotopt)
hold off

%% play with gui
app = networkGraphEdit('nt',NTg,'img',img)
