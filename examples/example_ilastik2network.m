% -----------
%% This is an example script to extract a network structure from an image
% The image should be grayscale: recommended to use a Probabilities file 
% obtained via Ilastik, where higher values correspond to pixels that
% are more likely to belong to the ER network
% When exporting from Ilastik, export only the label that corresponds to
% the ER
% ------------

%% Step 0:
% If you want the Ilastik step to go faster, first crop your movie to 
% a box around the cell. The below code can also do cropping in the
% post-Ilastik step just in case this wasn't done earlier

%% Step 1: Ilastik
% Open up Ilastik:
% 1) Load in raw data file (multipage tiff)
% 2) Feature Selection: select everything
% 3) Training: mark Label 1 as pixels belonging to ER and Label 2 as pixels
% not belonging to ER. Turn on Live Update to check on the prediction. 
% Make sure to train on several frames of the movie
% Turn on the eye next to Uncertainties to see which regions it is confused
% about.
% 4) Prediction export: Source = Probabilities. 
% Choose Export Image Settings, find the table that says Cutout Subregion and has rows labeled "x, y,
% z, c". On the "c" row, uncheck All, and select start=0, stop = 1. This
% will export only the Label 1 probabilities
% Format = multipage tiff
% 5) Remember to hit export button.

% ------------------
%% Step 2: Load in probabilities files and rawimage files into Matlab
% -------------------

% update dirname to be the directory where you are storing your files
dirname = '/data/proj/ERnetworks/Pekkurnaz/Alexa20250522/';
% update file name to be the probabilities file output by ilastik
probimgfile = [dirname 'COS7P13_mitoBFP_eGFP_mcherry-Sec61b_crop_Probabilities.tiff'];
% update file name to be the raw image file
rawimgfile = [dirname 'COS7P13_mitoBFP_eGFP_mcherry-Sec61b_crop.tiff']

info = imfinfo(rawimgfile);
pxperum = info(1).XResolution; % this is the resolution in pixels per micron

nframe = length(info); % total number of frames in movies
% if you want to work with only the first few frames (eg: for testing), uncomment the line
% below
nframe = 2;


probimgs = zeros(info(1).Height,info(1).Width,nframe);
rawimgs = probimgs;
for fc = 1:nframe
    disp([fc nframe])
    probimgs(:,:,fc) = imread(probimgfile,fc);
    rawimgs(:,:,fc) = imread(rawimgfile,fc);
end

%% threshold the first image and use this to define a cropping box
docrop = true; % do you want to do cropping?
% expand the cropping box by some number of pixels on either side
pxexpandbox = 30;

if (docrop)
    thresh = graythresh(probimgs(:,:,1));
    BW = imbinarize(probimgs(:,:,1),thresh);

    % define cropping box
    CC = bwconncomp(BW);
    props = regionprops(CC,"Area","BoundingBox");
    [maxArea,maxIdx] = max([props.Area]);

    cropbox = props(maxIdx).BoundingBox;

    cropbox = [max(cropbox(1:2)-pxexpandbox,0) min(cropbox(3:4)+2*pxexpandbox,fliplr(size(BW)))];

    % image of just the largest component
    BWcon = 0*BW;
    BWcon(CC.PixelIdxList{maxIdx}) = 1;

    figure
    imshow(BWcon,[])
    h_crop = drawrectangle('Position',cropbox)
    title('Thresholded Frame 1 and cropping box')
else
    docrop = [0 0 fliplr(size(probimgs(:,:,1)))];
end

%% Crop all raw images and probability images.
if docrop
    for fc = 1:nframe
        disp(fc)
        imgcrop = imcrop(probimgs(:,:,fc),cropbox);
        probimgcrop(:,:,fc) = imgcrop;
        rawimgcrop(:,:,fc) = imcrop(mat2gray(rawimgs(:,:,fc)),cropbox);
    end
else
    probimgcrop = probimgs;
    rawimgcrop = rawimgs;
end

%% Binarize (segment) the probability images
BWimgs = probimgcrop;
for fc = 1:nframe
    thresh = graythresh(probimgcrop(:,:,fc));
    BWimgs(:,:,fc) = imbinarize(probimgcrop(:,:,fc),thresh); 
end

%% Manually identify the perinuclear region that you do *not* want to work with
% Click on the figure to draw a polygon. It will stop drawing when you
% close the polygon by clicking on the original point.
imscl = imadjust(rawimgcrop(:,:,1), [0.05,0.4], [0,1], 0.4);
imshow(imscl,[])
% draw perinuclear region you don't want to extract
nuch = drawpolygon("Color","c")

%% Manually adjust the polygon contour if desired. Then extract its coordinates
nuccoords = nuch.Position;

%% black out perinuclear regions and drop all connected comments that are too small
% minimal area of connected component (in pixels)
% drop anything smaller than this
mincomponentarea = 500; 

nucmask = poly2mask(nuccoords(:,1),nuccoords(:,2),size(BWimgs,1),size(BWimgs,2));
nucmask = 1-nucmask;

for fc = 1:nframe
    img = nucmask.*BWimgs(:,:,fc);

    CC = bwconncomp(img);
    props = regionprops(CC);
    areas = [props.Area];
    dropcc = find(areas<mincomponentarea);
    droppx = vertcat(CC.PixelIdxList{dropcc}); % these pixels are zeroed out

    img(droppx) = 0;
    BWimgs(:,:,fc) = img;
end

% visualize 1 frame
imshow(BWimgs(:,:,1))

%% Extract networks from BW images
% Warning: This step is very slow! Takes about 1 min per frame.
clear allnetworks
options = struct('dodisplay',1);
for fc = 1:nframe
    disp(sprintf('Working on frame %d',fc));
    tic
    [NT,skelimage,opt] = getNetworkFromBWImage(BWimgs(:,:,fc),options);
    toc
    allnetworks(fc) = NT;
end

%% Save network objects for later use (to be reloaded into matlab sa needed)
savefile = [dirname 'COS7P13_mitoBFP_eGFP_mcherry-Sec61b_networks.mat'];
save(savefile,'dirname','probimgfile','rawimgfile','allnetworks','cropbox')

%% Output info on the networks to be loaded in other software (eg: with python, etc)
% one text file per network
% set the below flag to true if you want to output coordinates in terms of
% the original image file (as opposed to the cropped one)
% Output file format is:
% (1) lines that give you the n-th node x and y coordinates and some value v
% associated with the node (default is just 0)
% NODE n x y v 
% (2) lines that list the edges of the network. Edge E connects nodes n1
% and n2, and has length L
% EDGE E n1 n2 L
% (3) Lines that list coordinates of the path tracing along each edge. Edge
% E has N tracing points along it with coordinates x_1,x_2,...,x_N,
% y_1,y_2,..y_N
% EDGEPATH E N x_1 x_2 ... x_N y_1 y_2 ... y_N

shift_to_orig_coords = true;

outoptions = struct('WRITEPATHS',true);
for fc = 1:length(allnetworks)
    outfile = [dirname sprintf('COS7P13_mitoBFP_eGFP_mcherry-Sec61b_network_%d.txt',fc)];
    disp(sprintf('Outputting to file %s',outfile))
    NT = copy(allnetworks(fc));

    if (shift_to_orig_coords)
        % shift to original coordinates
        NT.nodepos = NT.nodepos + cropbox(1:2);
        for ec = 1:NT.nedge
            NT.edgepath{ec} = NT.edgepath{ec} + cropbox(1:2);
        end
    end

    NT.outputNetwork(outfile,outoptions);
end