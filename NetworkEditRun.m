function NetworkEditRun
    addpath('/home/matlab/Lena/networktools/');
    addpath('/home/matlab/Lena/networktools/gui/');
    addpath('/home/matlab/Lena/networktools/examples/');
    addpath('/home/matlab/Lena/');

    [file,path] = uigetfile('*.mat');
    load([path file]);
    networkEdit(NT, img, plotopt)
return

