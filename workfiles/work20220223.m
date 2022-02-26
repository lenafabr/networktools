% test gui
load('./examples/exampleERnetwork.mat')

%%
plotopt = struct('nodecolor',[1 0 0],'nodesize',20)
plotopt.edgeplotopt = {'LineWidth',2,'Color','g'};

networkEdit('NT',NT,'img',img,'plotopt',plotopt)