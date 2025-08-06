%% Load in probabilities file
dirname = '/data/proj/ERnetworks/Pekkurnaz/Alexa20250522/';
probimgfile = [dirname 'COS7P13_mitoBFP_eGFP_mcherry-Sec61b_Probabilities.tiff'];
rawimgfile = [dirname 'COS7P13_mitoBFP_eGFP_mcherry-Sec61b.tif']

info = imfinfo(probimgfile);
nframe = length(info);
%nframe = 10; % work with just a few frames to start with

probimgs = zeros(info(1).Width,info(1).Height,length(info));
rawimgs = probimgs;
for fc = 1:nframe
    disp([fc nframe])
    probimgs(:,:,fc) = imread(probimgfile,fc);
    rawimgs(:,:,fc) = imread(rawimgfile,fc);
end


%%
info = imfinfo(rawimgfile);
pxperum = info(1).XResolution;
%% threshold the first image
thresh = graythresh(probimgs(:,:,1));
BW = imbinarize(probimgs(:,:,1),thresh);
imshowpair(probimgs(:,:,1),BW)

%% define cropping box 
CC = bwconncomp(BW);
props = regionprops(CC,"Area","BoundingBox");
[maxArea,maxIdx] = max([props.Area]);

cropbox = props(maxIdx).BoundingBox;
% expand the cropping box by some number of pixels
pxexpandbox = 30;
cropbox = [cropbox(1:2)-pxexpandbox cropbox(3:4)+2*pxexpandbox];

% image of just the largest component
BWcon = 0*BW;
BWcon(CC.PixelIdxList{maxIdx}) = 1;

imshow(BWcon,[])
drawrectangle('Position',cropbox)

%% crop all images
clear BWimgs
for fc = 1:nframe
    disp(fc)
    imgcrop = imcrop(probimgs(:,:,fc),cropbox);
    thresh = graythresh(imgcrop);
    BWimgs(:,:,fc) = imbinarize(imgcrop,thresh); 
    rawimgcrop(:,:,fc) = imcrop(mat2gray(rawimgs(:,:,fc)),cropbox);
end
%% visualize
imscl = imadjust(rawimgcrop(:,:,1), [0.05,0.4], [0,1], 0.4);
imshowpair(imscl,BWimgs(:,:,1))

%% Save cropped raw images for later use
tiffile = [dirname 'COS7P13_mitoBFP_eGFP_mcherry-Sec61b_crop.tiff'];
imwrite(rawimgcrop(:,:,1),tiffile);
for fc = 2:nframe
    imwrite(rawimgcrop(:,:,fc),tiffile,"WriteMode","append");
end

%% Manually identify perinuclear region
imshow(imscl,[])
% draw perinuclear region you don't want to extract
nuch = drawpolygon("Color","c")
%%
nuccoords = nuch.Position;
%% black out perinuclear regions
nucmask = poly2mask(nuccoords(:,1),nuccoords(:,2),size(BWimgs,1),size(BWimgs,2));
nucmask = 1-nucmask;

for fc = 1:nframe
    BWimgs(:,:,fc) = nucmask.*BWimgs(:,:,fc);
end

%% Extract networks from BW images
clear allnetworks
options = struct('dodisplay',1);
for fc = 1:5%nframe
    disp(fc)
    [NT,skelimage,opt] = getNetworkFromBWImage(BWimgs(:,:,fc),options);
    allnetworks(fc) = NT;
end

%%
savefile = [dirname 'COS7P13_mitoBFP_eGFP_mcherry-Sec61b_networks.mat'];
save(savefile,'dirname','probimgfile','rawimgfile','allnetworks','cropbox')

%% Clean networks by removing very short tips
mintiplen = 0.5*pxperum; % minimal length of terminal edge to keep

for nc = 1:length(allnetworks);
    NT = allnetworks(nc);
    NT.mergeAllEdgePaths();

    tips = find(NT.degrees==1);
    elen = NT.edgelens(NT.nodeedges(tips,1));
    shorttips = find(elen<mintiplen);
    rmnodes = tips(shorttips);
    dokeep = true(1,NT.nnode);
    dokeep(rmnodes) = false;
    NT.keepNodes(find(dokeep));
end

%% remove break up any single-edge loops
for nc = 1:length(allnetworks)
    NT = allnetworks(nc);
    for ec = 1:NT.nedge
        if (NT.edgenodes(ec,1)==NT.edgenodes(ec,2))
            disp(sprintf('Breadking edge %d on network %d', ec, nc))
            NT.breakEdge(ec,0.5,true)
        end
    end    
end

%%
savefile = [dirname 'COS7P13_mitoBFP_eGFP_mcherry-Sec61b_networks_cleaned.mat'];
save(savefile,'dirname','probimgfile','rawimgfile','allnetworks','cropbox','mintiplen')

%% reload images and networks
load([dirname 'COS7P13_mitoBFP_eGFP_mcherry-Sec61b_networks_cleaned.mat']);
imgcropfile = [dirname 'COS7P13_mitoBFP_eGFP_mcherry-Sec61b_crop.tiff'];
clear rawimgcrop
for fc = 1:nframe
    tmp = imread(imgcropfile,fc);    
    % normalize
    rawimgcrop(:,:,fc) = double(tmp)/max(double(tmp(:)));
end

%% plot network over image

img = imadjust(rawimgcrop(:,:,1),[0.05,0.5],[0, 1],0.5);
imshow(img(:,:,1))
hold all
allnetworks(1).plotNetwork(struct('plotoverimage',true))
hold off

%%
addpath(genpath('~/proj/networktools/gui'))
%% adjust network if needed (old gui)
% plotting options (avoids drawing black lines on black background)
NT = allnetworks(1);
plotopt = struct('nodesize',20,'nodecolor',[1 0 0],'datatipindex',true)
plotopt.edgeplotopt = {'LineWidth',2,'Color','g'};

networkEdit('NT',NT,'img',img,'plotopt',plotopt)

%% adjust network (new gui)
NTg = NetworkGraphObj();
NTg.getFromOldNetwork(allnetworks(1));
NTg0 = copy(NTg);

%%
app = networkGraphEdit('nt',NTg,'img',img)


% -------------------
%% Automated tracking of tips (does not work well for growing tips
% -------------------
%% Find utrack code for tracking tips
utrackdir='~/proj/minNetworkDynamics/u-track-master/';
addpath(utrackdir)
addpath(genpath(utrackdir))

NT = allnetworks;

%% (2) create movieInfo structure
movieInfo=struct();
movieInfo.xCoord=[];
movieInfo.yCoord=[];
movieInfo.amp=[];
LCUTOFF=.1;


for fm=1:nframe
    fm    
    % find all d1nodes and their connected edge lengths
    d1nodes=find(NT(fm).degrees==1);
    d1edges=NT(fm).nodeedges(d1nodes);
    d1lengths=NT(fm).edgelens(d1edges);

    inds=find(d1lengths>LCUTOFF);

    ld1s=d1nodes(inds);
    lenld1s=NT(fm).edgelens(ld1s);
    xs = NT(fm).nodepos(ld1s,1);
    ys = NT(fm).nodepos(ld1s,2);
    
    % set movieInfo, x, y coordinates and amplitudes, setting all amps to 1
    movieInfo(fm).xCoord = [xs zeros(size(xs))]; 
    movieInfo(fm).yCoord = [ys zeros(size(ys))];
    movieInfo(fm).amp = [ones(size(lenld1s)) zeros(size(lenld1s))];
end

%% (3) u-track
run scriptTrack1DNodes.m

%% (4) find long tracks, within eroded boundary of cell (since lots of d1 nodes appear at boundary

% parameter to play with!
bndErode = -.1;

x=NT(fm).nodepos(:,1);y=NT(fm).nodepos(:,2);
k=boundary(x,y);
pbnd=polyshape(x(k),y(k));
pbnd2=polybuffer(pbnd,bndErode);

TL=[];
for tc=1:length(tracksFinal)
    TL(tc) = length(tracksFinal(tc).tracksFeatIndxCG);
    tC=tracksFinal(tc).tracksCoordAmpCG;
    st(tc) = isinterior(pbnd2,tC(1),tC(2));
    en(tc) = isinterior(pbnd2,tC(end-7),tC(end-6));
end
inds1=find((TL>=3));
inds2=find(st);
inds3=find(en);
inds=intersect(inds1,inds2);
inds=intersect(inds,inds3);

% final set of tracks
tracks=tracksFinal(inds);

%% (5) visualize
options=struct();
options.ERdata = [dirname 'COS7P13_mitoBFP_eGFP_mcherry-Sec61b_crop.tiff'];
options.frames = 1:10;
options.axlim=[0 cropbox(3) 0 cropbox(4)];
options.overlay=true;
options.getMovie=false;
img = imread(options.ERdata,1);
options.clim = [0,max(img(:))*0.3];
% options.frames=1:10;
% options.dir='/home/zscott/Desktop/';
plotGrowthTracks(tracks,options);

%%
savefile = [dirname 'COS7P13_mitoBFP_eGFP_mcherry-Sec61b_networks_tracks.mat'];
save(savefile,'dirname','probimgfile','rawimgfile','allnetworks','cropbox','tracks','options')


% -----------
%% transition to manual track curation
% -----------


%% load in images
imgfile = [dirname 'COS7P13_mitoBFP_eGFP_mcherry-Sec61b_crop.tiff'];
imginfo = imfinfo(imgfile);
nframe = length(imginfo); % number of frames

% load all the frames
for sc = 1:nframe
    imgs(:,:,sc) = imread(imgfile,sc);  
    imgs(:,:,sc) = imadjust(imgs(:,:,sc),[0.05,0.4],[0, 1],0.4);
end


%% set up data structures for tracks

pxperum = 1%imginfo(1).XResolution;

[tracks,framepts] = setTrackData(tracks,nframe,pxperum);
ntrack = length(tracks);
nframe = length(framepts);

%% check out 1 track
allnetworks(1).plotNetwork()
set(gca,'YDir','reverse')

hold all
plot(tracks(1).pos(:,1),tracks(1).pos(:,2),'r.-')
hold off

%% Run GUI
app =tracktubules(imgs,framepts,allnetworks);
app.cleanUpPositionBug()
%% get new tracks
newtracks = framePts2Tracks(app.framepts,app.UITable.Data,tracks,pxperum)

%% check that everything is consistent
[tracks2,framepts2] = setTrackData(newtracks,nframe,pxperum)

app2 =tracktubules(imgs,framepts2);
app2.cleanUpPositionBug();


%% Invert probabilities and crop
probimgscrop = zeros(size(rawimgcrop));
for fc = 1:nframe
    probimgscrop(:,:,fc) = imcrop(probimgs(:,:,fc),cropbox);
    probimgscrop(:,:,fc) = nucmask.*(1-probimgscrop(:,:,fc));
end

%% Save inverted probabilities
probtiffile = [dirname 'COS7P13_mitoBFP_eGFP_mcherry-Sec61b_probinv_crop.tiff'];
imwrite(probimgscrop(:,:,1),probtiffile);
for fc = 2:nframe
    imwrite(probimgscrop(:,:,fc),probtiffile,"WriteMode","append");
end

%% reload some info
dirname = '/data/proj/ERnetworks/Pekkurnaz/Alexa20250522/';
rawimgfile = [dirname 'COS7P13_mitoBFP_eGFP_mcherry-Sec61b.tif']
info = imfinfo(rawimgfile);
pxperum = info(1).XResolution;

%% Load in data on polygon objects
data = readtable([dirname 'COS7P13_mitoBFP_eGFP_mcherry-Sec61b_crop_table.csv']);
%%
% pick objects more likely good than bad;
thresh = 0.5;
goodind = find(data.ProbabilityOfLabel1>thresh);
polyareas = data.SizeInPixels(goodind)/pxperum^2;
meanA = mean(polyareas);

histogram(log10(polyareas/meanA),100,'Normalization','pdf')
set(gca,'YScale','log')
hold all
xvals = logspace(-2,0);
%semilogy(log10(xvals),200*sqrt(xvals.^2),'r')
%hold off
xlim([-2,1])

%% Zuben's statistics for experimental data
load('../../results/areaDistsForGreg.mat')

%%
histogram(log10(allAreasC/mean(allAreasC)),'Normalization','pdf')
set(gca,'YScale','log')
hold all
xvals = logspace(-2,0);
semilogy(log10(xvals),sqrt(xvals),'r')
hold off
xlim([-2,1])

%% crop to desired small box
for fc = 1:size(rawimgs,3)
    imgsq(:,:,fc) = imcrop(rawimgs(:,:,fc),h.Position);
end
%%
app =tracktubules(imgsq,framepts,allnetworks);
app.cleanUpPositionBug()