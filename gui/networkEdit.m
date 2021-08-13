function varargout = LenaFig(varargin)
% LENAFIG MATLAB code for LenaFig.fig
%      LENAFIG, by itself, creates a new LENAFIG or raises the existing
%      singleton*.
%
%      H = LENAFIG returns the handle to a new LENAFIG or the handle to
%      the existing singleton*.
%
%      LENAFIG('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in LENAFIG.M with the given input arguments.
%
%      LENAFIG('Property','Value',...) creates a new LENAFIG or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before LenaFig_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to LenaFig_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help LenaFig

% Last Modified by GUIDE v2.5 12-Aug-2021 12:54:58

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @LenaFig_OpeningFcn, ...
                   'gui_OutputFcn',  @LenaFig_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before LenaFig is made visible.
function LenaFig_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to LenaFig (see VARARGIN)

% Choose default command line output for LenaFig
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
cd /home/matlab//Lena/

% UIWAIT makes LenaFig wait for user response (see UIRESUME)
% uiwait(handles.mainFig);


% --- Outputs from this function are returned to the command line.
function varargout = LenaFig_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Manage menu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function menuLoad_Callback(hObject, eventdata, handles)
global newf NTobj imgObj plotoptObj selNodes selEdges 
    close(figure(1))
    newf =[];
    NTobj = [];
    imgObj = [];
    plotoptObj = [];
    selNodes = [];
    selEdges = [];

    load('./examples/exampleERnetwork')
    
    NTobj = NT;
    imgObj = img;
    plotoptObj = plotopt;
    
    dispNetWithImage();
return

function dispNetWithImage()
global newf NTobj imgObj plotoptObj nodeplotH edgeplotH imageH selNodes;
    try
        %close(newf);
        delete(nodeplotH)
        for lc = 1:length(edgeplotH)
            delete(edgeplotH(lc))
        end
    end

    newf=figure(1);
    
    if (exist('imageH','var'))
        if isempty(imageH)
            imageexists = false;
        else
            imageexists = isvalid(imageH);
        end
    else
        imageexists = false;
    end
    
    if (~imageexists) % redraw figure
        imageH = imshow(imgObj,[0,0.9]);
        set(newf, 'Position', [20 20 500 500]);
        set(gca,'Position',[0,0,1,1])
    end
    hold all

    % plot network
    newf=figure(1);
    hold on   
    %newf=figure(1);
    plotoptObj.datatipindex = true;
    [nodeplotH,edgeplotH] = NTobj.plotNetwork(plotoptObj);
    hold off
    
    if ~isempty(selNodes)
        figure(newf);
        hold on
        selscatH = scatter(NTobj.nodepos(selNodes,1), NTobj.nodepos(selNodes,2),...
            20, 'b', 'filled');
        selscatH.PickableParts='none';
        hold off       
    end

    nodeplotH.PickableParts = 'none';
    for lc = 1:length(edgeplotH); edgeplotH(lc).PickableParts = 'none'; end
    
%     % set data tip properties   
%     dt = nodeplotH.DataTipTemplate;
%     dt.DataTipRows(1).Value = 1:NTobj.nnode;
%     dt.DataTipRows(1).Label = '';
%     dt.DataTipRows(2:end) = [];
%     dt.FontSize=6;     
    if (~imageexists)
        set(gca,'Position',[0,0,1,1])
    end
return

function menuClear_Callback(hObject, eventdata, handles)
    global newf NTobj selNodes 
    
    close(figure(1))
    newf =[];
    NTobj = [];
    selNodes = [];
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Manage nodes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ind = findNearestNode(nodepos, xy)
    dNodes = nodepos - xy;
    d2 = dNodes(:,1).^2 + dNodes(:,2).^2;
    M = min(d2);
    ind = find(d2 == M);       
return

function ind = selectNode(addSelected, color)
    global newf NTobj nodeplotH selNodes

    ind = [];
    nodeplotH.PickableParts = 'all';

    figure(newf)
    w = 0;
    while ~w
        w = waitforbuttonpress;
    end
    %input('On figure, select nodes for removal. Then press enter.')
    
    figure(newf);
    datatips = findobj(gca,'Type','datatip');
    
    if ~isempty(datatips)
        ind = [datatips.DataIndex];
        scatter = findobj(gca,'Type','scatter');
        for i=1:length(ind)
            scatter.CData(ind(i),:) = color;
        end
        
        if addSelected
            selNodes = [selNodes ind];
            selNodes = unique(selNodes);
        end
        
        delete(datatips)
    end
    nodeplotH.PickableParts = 'none';
return
    
function pushbuttonSelectNode_Callback(hObject, eventdata, handles)  
    ind = selectNode(true, [0 0 1]);
return

function pushbuttonUnselectNode_Callback(hObject, eventdata, handles)
    global selNodes

    iSel = selectNode(false, [1 0 0]);
    
    if ~isempty(iSel)
         [val, iRem, iInd] = intersect(selNodes, iSel);
         selNodes(iRem) = [];         
    end
return

function pushbuttonUnselectAllNodes_Callback(hObject, eventdata, handles)
    global newf selNodes
        
    figure(newf)
    scatter = findobj(gca,'Type','scatter');
    for i=1:length(selNodes)
        scatter.CData(selNodes(i),:) = [1 0 0];
    end
    selNodes = [];
    
return

function pushbuttonSelArea_Callback(hObject, eventdata, handles)
    global newf NTobj selNodes selEdges
    selEdges = unique(selEdges);
    figure(newf);
    hold on
    
    roi = drawpolyline;
    xp = roi.Position(:,1);
    yp = roi.Position(:,2);
    xc = [xp' roi.Position(1,1)]';
    yc = [yp' roi.Position(1,2)]';
    plot(xp,yp, 'm')
    
    x = NTobj.nodepos(:,1);
    y = NTobj.nodepos(:,2);
    [in on] = inpolygon(x,y, xp,yp);
    delete(roi);
    hold off
    
    ind = find(in);
    scatter = findobj(gca,'Type','scatter');
    for i=1:length(ind)
        scatter.CData(ind(i),:) = [0 0 1];
    end

    selNodes = [selNodes' ind']';
    selNodes = unique(selNodes);
    
    %removeSelected();
return

function pushbuttonAddNode_Callback(hObject, eventdata, handles)
    global newf NTobj selNodes
    
    figure(newf);
    roi = drawpoint;
    xy = roi.Position;
    delete(roi);
    
    NTobj.nnode = NTobj.nnode + 1;
    nnode = NTobj.nnode;
    NTobj.nodepos(nnode,:) = xy';
    NTobj.degrees(nnode) = 0;
    NTobj.nodenodes(nnode,:) = [0 0 0 0];
    NTobj.nodeedges(nnode,:) = [0 0 0 0];
        
    dispNetWithImage();
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Manage edges
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function iSel = selectEdge(color)
    global newf selEdges edgeplotH
              
    iSel=[];
    for ie = 1:length(edgeplotH); 
        edgeplotH(ie).PickableParts = 'all'; 
    end
   
    figure(newf)
    w = 0;
    while ~w
        w = waitforbuttonpress;
    end
    %input('On figure, select edges. Then press enter.')
    
    figure(newf)
    datatips = findobj(gca,'Type','datatip');
    
    if ~isempty(datatips)
        for idt=1:length(datatips)
            Ln = datatips(idt).Parent;
            Ln.Color = color;
            iSel = [iSel Ln.edgeind];
        end
        
        delete(datatips)
        selEdges = [selEdges iSel];
    end

    for lc = 1:length(edgeplotH); 
        edgeplotH(lc).PickableParts = 'none'; 
    end
return

function pushbuttonSelectEdge_Callback(hObject, eventdata, handles)
    iSel = selectEdge('b');
return

function pushbuttonUnselectEdge_Callback(hObject, eventdata, handles)
    global selEdges
    
    iSel = selectEdge('g');
    
    if ~isempty(iSel)
        [val, iRem, iInd] = intersect(selEdges, iSel);
        selEdges(iRem) = [];
    end
return

function pushbuttonUnselectAllEdges_Callback(hObject, eventdata, handles)
    global newf selEdges

    figure(newf)
    L = findobj(gca,'Type','line');
    for i=1:length(L)
        L(i).Color = 'g';
    end
    selEdges = [];
return

function pushbuttonAddEdge_Callback(hObject, eventdata, handles)
    global newf NTobj edgeplotH
    
    figure(newf);
    
    roi = drawpolyline;
    pnt = roi.Position;
    delete(roi);
    
    ind1 = findNearestNode(NTobj.nodepos, pnt(1,:));
    pnt(1,:) = NTobj.nodepos(ind1,:);
    ind2 = findNearestNode(NTobj.nodepos, pnt(end,:));        
    pnt(end,:) = NTobj.nodepos(ind2,:);
    x = pnt(:,1);
    y = pnt(:,2);
    
    i = NTobj.nedge + 1;
    NTobj.nedge = i;
    
    NTobj.edgenodes(i,:) = [ind1 ind2];
    NTobj.edgepath{i} = [x y];
    NTobj.cumedgelen{i} =[];
    
    %setCumEdgeLen(NTobj, i,true);
    
    CLF = hypot(diff(x), diff(y));   % Calculate integrand from x,y derivatives
    NTobj.edgelens(i) = trapz(CLF); % Integrate to calculate arc length
    
    NTobj.edgeedges(i,:,:) = NTobj.edgeedges(i-1,:,:);
    
    hold on
    edgeplotH(i) = plot(x,y, 'Color','g', 'LineWidth',2);
    hold off
    edgeplotH(i).addprop('edgeind');
    edgeplotH(i).edgeind = i;
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function removeSelected_unfinished()
global NTobj newf selNodes selEdges

if isempty(selNodes) & isempty(selEdges)
    return
end

dokeep = true(1,NTobj.nnode);
dokeep(selNodes) = false;
keepind = find(dokeep);
keepNodes(NTobj,keepind);


selNodes=[];

%     if ~isempty(selEdges)
%         NTobj.edgepath(selEdges) =[];
%         NTobj.cumedgelen(selEdges) =[];
%         NTobj.edgelens(selEdges) =[];
%         NTobj.edgeedges(selEdges,:,:) =[];
%         NTobj.edgenodes(selEdges,:) =[];
%
%         NTobj.nedge = NTobj.nedge - length(selEdges);
%         selEdges = [];
%     end

redraw();
return

function removeSelected()
    global NTobj newf selNodes selEdges
    
    if isempty(selNodes) & isempty(selEdges)
        return
    end
    
    figure(newf)
    if ~isempty(selNodes)
        for i=1:length(selNodes)
            iSel = selNodes(i);
            
            % Remove edges stemming from this node
            for j=1:NTobj.nnode
                iMatch = find(NTobj.nodenodes(j,:)==iSel);
                if ~isempty(iMatch)
                    buf = NTobj.nodenodes(j,:);
                    buf(iMatch) = [];
                    NTobj.nodenodes(j,:) = [buf 0];
                end
            end
            
            for j=1:NTobj.nedge
                iMatch = find(NTobj.edgenodes(j,:)==iSel);
                if ~isempty(iMatch)
                    selEdges = [selEdges j];
                end
            end
        end
        selEdges = sort(unique(selEdges), 'descend');
        
        NTobj.degrees(selNodes) = [];
        NTobj.nodepos(selNodes,:) = [];
        NTobj.nodenodes(selNodes,:);
        NTobj.nodeedges(selNodes,:) = [];
        
        NTobj.nnode = NTobj.nnode - length(selNodes);
        selNodes=[];
    end
    
    if ~isempty(selEdges)
        NTobj.edgepath(selEdges) =[];
        NTobj.cumedgelen(selEdges) =[];
        NTobj.edgelens(selEdges) =[];
        NTobj.edgeedges(selEdges,:,:) =[];
        NTobj.edgenodes(selEdges,:) =[];
        
        NTobj.nedge = NTobj.nedge - length(selEdges);
        selEdges = [];
    end
    
    redraw();
return

function pushbuttonRemoveSelected_Callback(hObject, eventdata, handles)
    removeSelected();
return

function redraw()
    global newf
    figure(newf);
    
    scatter = findobj(gca,'Type','scatter');
    delete(scatter);
    
    Lines = findobj(gca,'Type','line');
    delete(Lines);
    
    dispNetWithImage();
return