load('./workfiles/example_network_Halo_Sec61_4b.mat')

%% plot image
imshow(img,[0,0.9])
hold all

% plot network
NT.plotNetwork(plotopt)
hold off