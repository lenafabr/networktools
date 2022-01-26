%% testing network guis

load('/UCSD/proj/dendriticmito/results/workspace_filteredNT_20210817.mat')
imgfilename = '~/proj/dendriticmito/data/MCFO-HSN-A61-2.tif';
img = imread(imgfilename);
filteredNT.edgevals={};
for ec = 1:filteredNT.nedge
    filteredNT.edgevals{ec} = [];
end

plotopt = struct('nodecolor',[1 0 0],'nodesize',20)
plotopt.edgeplotopt = {'LineWidth',2,'Color','g'};

%% make new network

%%
networkEdit('NT',filteredNT,'img',img,'plotopt',plotopt)