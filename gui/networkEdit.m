function varargout = networkEdit(varargin)
% networkEdit MATLAB code for networkEdit.fig
%      networkEdit, by itself, creates a new networkEdit or raises the existing
%      singleton*.
%
%      H = networkEdit returns the handle to a new networkEdit or the handle to
%      the existing singleton*.
%
%      networkEdit('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in networkEdit.M with the given input arguments.
%
%      networkEdit('Property','Value',...) creates a new networkEdit or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before networkEdit_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to networkEdit_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help networkEdit

% Last Modified by GUIDE v2.5 19-Sep-2023 15:21:56

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @networkEdit_OpeningFcn, ...
                   'gui_OutputFcn',  @networkEdit_OutputFcn, ...
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


% --- Executes just before networkEdit is made visible.
function networkEdit_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to networkEdit (see VARARGIN)
    global newf NTobj NTlocal imgObj plotoptObj selNodes selEdges guilock
    % Choose default command line output for networkEdit
    handles.output = hObject;
    
    % Update handles structure
    guidata(hObject, handles);

    % specific for AF
%     addpath /home/matlab/Lena/networktools;
%     addpath /home/matlab/Lena/networktools/gui;
    
%% set defaults
    NTobj = [];
    imgObj = [];
    plotoptObj = struct();

    for index = 1:2:length(varargin)
        switch(lower(varargin{index}))
            case 'nt'
                NTobj = varargin{index+1};
            case 'img'
                imgObj = varargin{index+1};
            case('plotopt')
                plotoptObj = varargin{index+1};
        end        
    end
    if ~isempty(NTobj)
        NTlocal = getNTlocal(NTobj);
    end
    
    newf =[];
    selNodes = [];
    selEdges = [];
    if (~isempty(NTobj))
        if (isempty(NTobj.edgewidth)) 
            NTobj.edgewidth = cell(NTobj.nedge,1);
        end
        % reset edgevals
        NTobj.edgevals={};
        for ec = 1:NTobj.nedge
            NTobj.edgevals{ec} = [];
        end
        dispNetWithImage();
    end

    hSlider1 = findobj('Tag', 'sliderContrast');
    hSlider1.Min = 0;
    hSlider1.Max = 4;
    hSlider2 = findobj('Tag', 'sliderBrightness');
    hSlider2.Min = -0.5;
    hSlider2.Max = 0.5;
    guilock = false;

    if (isempty(NTobj))
        ActionsDisable(handles);
    else
        ActionsEnable(handles);
    end
    
    set(handles.pushbuttonWidthCalc,'Enable','off');
    set(handles.pushbuttonWidthCancel,'Enable','off');
    % UIWAIT makes networkEdit wait for user response (see UIRESUME)
% uiwait(handles.mainFig);

return
    
% --- Outputs from this function are returned to the command line.
function varargout = networkEdit_OutputFcn(hObject, eventdata, handles) 
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
function LoadData()
    global newf NTobj NTlocal imgObj plotoptObj selNodes selEdges fileName imageH
%     global imgCData0
    
    close(figure(1))
    newf =[];
    NTobj = [];
    imgObj = [];
    plotoptObj = [];
    selNodes = [];
    selEdges = [];
    
    [file,path] = uigetfile('../*.mat');
    fileName = [path file];
% fileName = '/home/matlab/Lena/networktools/exampleERnetwork.mat';
    load(fileName);

    NTobj = NT;
    imgObj = img;
    plotoptObj = plotopt;
    if isempty(NTobj.edgewidth)
        NTobj.edgewidth = cell(NT.nedge,1);
    end
    NTlocal = getNTlocal(NTobj);
    
    figure(1);       
    imageH = imshow(imgObj,[]);
%     imgCData0 = imageH.CData;

    dispNetWithImage();
    figure(newf);
    
return

function NTlocal = getNTlocal(NTobj)
    NTlocal.nodepos = NTobj.nodepos;
    NTlocal.edgenodes = NTobj.edgenodes;
    NTlocal.edgepath = NTobj.edgepath;
    NTlocal.nodeactive = true(NTobj.nnode,1);
    NTlocal.edgeactive = true(NTobj.nedge,1);
    NTlocal.edgewidth=NTobj.edgewidth;
    NTlocal.filtered = true;
return

function menuLoad_Callback(hObject, eventdata, handles)
    LoadData();
    ActionsEnable(handles);
return

function menuSave_Callback(hObject, eventdata, handles)
    global NTobj NTlocal imgObj plotoptObj fileName
    
    if ~NTlocal.filtered
        filterActiveNetwork(NTlocal.nodepos, NTlocal.edgenodes,...
        NTlocal.nodeactive, NTlocal.edgeactive, NTlocal.edgepath,NTobj);
        NTlocal.filtered = true;
    end
    
    NT = NTobj;
    img = imgObj;
    plotopt = plotoptObj;
    newFile = [fileName(1:end-4) '_modified.mat'];
    uisave({'NT','img','plotopt'}, newFile);
return

function menuClear_Callback(hObject, eventdata, handles)
    global newf NTobj NTlocal selNodes 
    
    close(figure(1))
    newf =[];
    NTobj = [];
    NTlocal = [];
    selNodes = [];
    
    hSlider1 = findobj('Tag', 'sliderContrast');
    hSlider1.Value = 1;
    hSlider2 = findobj('Tag', 'sliderBrightness');
    hSlider2.Value = 0;
    
    ActionsDisable(handles)
return

function menuQuit_Callback(hObject, eventdata, handles)
    close all
    clear all
return
step

function dispNetWithImage()
global newf NTobj imgObj plotoptObj nodeplotH edgeplotH imageH selNodes;
        
    try
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
            imageexists = isgraphics(imageH);
            if (imageexists)
                imageexists = isvalid(imageH);
            end
        end
    else
        imageexists = false;
    end
    
    if (~imageexists) % redraw figure
        if (~isempty(imgObj))
            % rescale image to be between 0 1 
            imgObj = double(imgObj)/max(double(imgObj(:))); 
            imageH = imshow(imgObj,[0 1]);
        end
        set(newf, 'Position', [40 40 500 500]);
        set(gca,'Position',[0,0,1,1])
    end
    hold all

    % plot network
    newf=figure(1);
    hold on   
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
    for lc = 1:length(edgeplotH); 
        edgeplotH(lc).PickableParts = 'none'; 
    end
    
    if (~imageexists)
        set(gca,'Position',[0,0,1,1])
    end
return

function sliderContrast_Callback(hObject, eventdata, handles)
    global imageH  imgObj guilock
    try
        a = get(hObject,'Value');
        hSlider2 = findobj('Tag', 'sliderBrightness');
        b = hSlider2.Value;
        imageH.CData = imgObj * a + b;
              
        imageH.CData(imageH.CData>1) = 1;
        imageH.CData(imageH.CData<0) = 0;
    catch exception
        guilock = false;
        disp(getReport(exception))
    end
return

function sliderBrightness_Callback(hObject, eventdata, handles)
    global imageH B imgObj guilock
    try
        hSlider1 = findobj('Tag', 'sliderContrast');
        a = hSlider1.Value;
        b = get(hObject,'Value');
        imageH.CData = imgObj*a+b
        
        imageH.CData(imageH.CData>1) = 1;
        imageH.CData(imageH.CData<0) = 0;
    catch exception
        guilock = false;
        disp(getReport(exception))
    end
return

function pushbuttonImgReset_Callback(hObject, eventdata, handles)
    global imageH imgObj

    imageH.CData = imgObj;
    hSlider1 = findobj('Tag', 'sliderContrast');
    hSlider1.Value = 1;
    hSlider2 = findobj('Tag', 'sliderBrightness');
    hSlider2.Value = 0;    
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Manage nodes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ind = selectNode(addSelected, color)
    global newf nodeplotH selNodes guilock

    disp('Use datatips to select desired nodes. Then hit any keyboard key (while the figure window is active).')
   
    try
        figure(newf)

        scatter = findobj(gca,'Type','scatter');
        PickableOn(scatter);
        PickableOn(nodeplotH);
        
        ind = [];
        w = 0;
        if ~exist('newf','var')
            return;
        end
        figure(newf);
        while ~w
            w = waitforbuttonpress;
        end
        figure(newf)
        datatips = findobj(gca,'Type','datatip');
        
        % Cancel selection if 'Esc' hit
        value = double(get(gcf,'CurrentCharacter'));
        if value==27
            figure(newf)
            delete(datatips)
            scatter = findobj(gca,'Type','scatter');
            PickableOff(nodeplotH);
            PickableOff(scatter);
            return
        end
                
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
    catch exception
        guilock = false;
        disp(getReport(exception))
    end
    if ~exist('newf','var')
        return;
    end
    figure(newf)
    scatter = findobj(gca,'Type','scatter');
    PickableOff(scatter);
    PickableOff(nodeplotH);
return
    
function pushbuttonSelectNode_Callback(hObject, eventdata, handles)
    global guilock
    if (guilock)
        disp('Cannot select nodes, gui is locked. Finish previous operation.')
        return
    end
    StartAction(handles, "Use datatips to select nodes. Then hit Enter to complete, Esc to cancel")
    set(gcf,'Pointer','Arrow');
    ind = selectNode(true, [0 0 1]);
    EndAction(handles)
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
    global newf NTlocal selNodes selEdges
    
    selEdges = unique(selEdges);
    figure(newf);
    hold on
    
    roi = drawpolyline;
    xp = roi.Position(:,1);
    yp = roi.Position(:,2);
    xc = [xp' roi.Position(1,1)]';
    yc = [yp' roi.Position(1,2)]';
    %plot(xp,yp, 'm')
    
    x = NTlocal.nodepos(:,1);
    y = NTlocal.nodepos(:,2);
    [in on] = inpolygon(x,y, xp,yp);
    delete(roi);
    
    ind = find(in);
    scatter = findobj(gca,'Type','scatter');
    for i=1:length(ind)
        scatter.CData(ind(i),:) = [0 0 1];
    end

    selNodes = [selNodes ind'];
    selNodes = unique(selNodes);
    
    %removeSelected();
    hold off
return

function pushbuttonAddNode_Callback(hObject, eventdata, handles)
    global NTlocal guilock
    
    StartAction(handles, ...
        'Click on positions of desired nodes. To finish hit Esc.')
    try
        P = []; H=[];
        i=0;
        
        set(gcf,'CurrentCharacter',char(1));
        value = 0;
        while true
            h = drawpoint;                     
           
            value = double(get(gcf,'CurrentCharacter'));
            if value==27
                break;
            end
            
            i = i+1;
            P{i} = h.Position;
            H = [H h];
         end
       
        scatter = findobj(gca,'Type','scatter');
        for i=1:length(P)
            nnode = size(NTlocal.nodepos,1) + 1;
            NTlocal.nodepos(nnode,:) = P{i}';
            NTlocal.nodeactive(nnode) = true;

            scatter.SizeData = [scatter.SizeData scatter.SizeData(1)];
            scatter.XData = [scatter.XData P{i}(1)];
            scatter.YData = [scatter.YData P{i}(2)];
            scatter.CData = [scatter.CData' [1 0 0]']';
            
            scatter.DataTipTemplate.DataTipRows.Value =...
                [scatter.DataTipTemplate.DataTipRows.Value  ...
                (max(scatter.DataTipTemplate.DataTipRows.Value)+1)];

            P{i} = [];
            delete(H(i));
        end
        P=[]; H=[];
        
        L = findobj(gca,'Type','line');
        for i=1:length(L)
            L(i).PickableParts = 'none';
            L(i).HitTest = 'off';
        end
        
        S = findobj(gca,'Type','scatter');
        for i=1:length(S)
            S(i).PickableParts = 'none';
            S(i).HitTest = 'off';
        end
    catch exception
        guilock = false;
        disp(getReport(exception))
    end
    NTlocal.filtered = false;
    EndAction(handles)
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Manage edges
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function iSel = selectEdge(color)
    global newf edgeplotH guilock NTlocal
    
    figure(newf)
    L = findobj(gca,'Type','line');
    PickableOn(L);
    PickableOn(edgeplotH);
    try
        figure(newf)

        iSel=[];  
        w = 0;
        display('Select desired edges. Then press any keyboard key (while the figure window is active)')
        
        while ~w
            w = waitforbuttonpress;
        end

        figure(newf)
        datatips = findobj(gca,'Type','datatip');
        
        % Cancel selection if 'Esc' hit
        value = double(get(gcf,'CurrentCharacter'));
        if value==27
            delete(datatips)
            PickableOff(L);
            PickableOff(edgeplotH);
            return
        end

        if ~isempty(datatips)
            for idt=1:length(datatips)
                Ln = datatips(idt).Parent;
                Ln.Color = color;
                iSel = [iSel Ln.edgeind];
            end

            delete(datatips)
            %selEdges = [selEdges iSel];
        end
    catch exception
        guilock = false;
        disp(getReport(exception))
    end
    figure(newf)
    L = findobj(gca,'Type','line');
    PickableOff(L);
    PickableOff(edgeplotH);
    
%     Lines = findobj(gca,'Type','line');
return

function pushbuttonSelectEdge_Callback(hObject, eventdata, handles)
    global guilock selEdges NTlocal
    if (guilock)
        disp('Cannot select edges, gui is locked. Finish previous operation.')
        return
    end
    
    StartAction(handles, "Use datatips to select edges. Then hit Enter to complete, Esc to cancel")
    set(gcf,'Pointer','Arrow');
    iSel = selectEdge('b');
    selEdges = [selEdges iSel];
    selEdges = sort(unique(selEdges));
    EndAction(handles)
    
%     edge = selEdges(1);
%     Lines = findobj(gca,'Type','line');
%     lineInd = edge2lineInd(edge, Lines, NTlocal.edgepath);
    
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
    global newf NTlocal edgeplotH guilock
    
    if (guilock)
        disp('Cannot draw new edge, gui is locked. Finish previous operation')
        return
    end
    guilock = true;
    
    try
        figure(newf);

        roi = drawpolyline;

        w = 0;
        display('Draw new edge path. When done, do right mouse click. Then adjust path and hit any key to continue')
        while ~w
            w = waitforbuttonpress;
        end
        figure(newf);
        pnt = roi.Position;
        delete(roi);

        ind1 = findNearestNode(NTlocal.nodepos, pnt(1,:));
        pnt(1,:) = NTlocal.nodepos(ind1,:);
        ind2 = findNearestNode(NTlocal.nodepos, pnt(end,:));        
        pnt(end,:) = NTlocal.nodepos(ind2,:);
        x = pnt(:,1); y = pnt(:,2);
        
        i = size(NTlocal.edgenodes,1) + 1;

        NTlocal.edgenodes(i,:) = [ind1 ind2];
        NTlocal.edgepath{i} = [x y];  
        NTlocal.edgewidth{i} = [];
        NTlocal.edgeactive(i) = true;
        
        hold on
        edgeplotH(i) = plot(x,y, 'Color','g', 'LineWidth',2);
        hold off
        edgeplotH(i).addprop('edgeind');
        edgeplotH(i).edgeind = i;
    catch exception
        guilock = false;
        disp(getReport(exception))
    end
    NTlocal.filtered = false;
    guilock = false;
return
 
function ind = findNearestNode(nodepos, xy)
    global NTlocal
    dNodes = nodepos - xy;
    d2 = dNodes(:,1).^2 + dNodes(:,2).^2;
    M = min(d2(NTlocal.nodeactive));
    ind = find(d2 == M);       
return

function iSel = findNearestEdge(xy)
    global NTobj NTlocal

    dMin = 1E20;
    iSel =-1;
    for i=1:NTobj.nedge
        if ~NTlocal.edgeactive
            continue;
        end
        d = NTobj.edgepath{i} - xy;
        d2 = d(:,1).^2 + d(:,2).^2;
        m = min(d2);
        if m<dMin
            iSel = i;
            dMin = m;
        end
    end
return

function [w, d, iSel] = getWidth(ep)
    global NTobj
    
    w = 0;
    d = 0;
    iSel = -1;

    X1 = ep(1,1); Y1 = ep(1,2); X2 = ep(2,1); Y2 = ep(2,2);
    xy = [0.5*(X1+X2) 0.5*(Y1+Y2)];
    % find nearest edge to center point of drawn segment
    iSel = findNearestEdge(xy);
    edge = NTobj.edgepath{iSel}; % path of the edge
    cum = NTobj.cumedgelen{iSel};
    
    % find intersection between drawn segment and edge
    seg1 = ep;
    segdiff = seg1(2,:)-seg1(1,:);
    w = norm(segdiff);
    
    for j = 1:size(edge,1)-1
        seg2 = edge(j:j+1,:);
        [t1,t2,intpt,doint] = segsegintersect(seg1', seg2');
        
        if doint
            xCros = intpt(1); yCros = intpt(2);
            d = cum(j)+norm(intpt'-edge(j,:));
            break
        end
    end
return

function pushbuttonWidthMeas_Callback(hObject, eventdata, handles)
    global newf NTobj guilock Hlist
              
    if (guilock)
        disp('Cannot measure widths, gui is locked. Finish previous operation.')
        ind = [];
        return
    end   
    StartAction(handles, 'Measure width of desired edges. Press Esc when done.')
    %Calculating edge width...')

    Hlist = []; 
    % Restore saved interactive lines
    figure(newf)
    for i=1:length(NTobj.edgewidth)
        if isempty(NTobj.edgewidth{i})
            continue;
        end
        for j=1:size(NTobj.edgewidth{i}, 1)
            L = images.roi.Line(gca);
            X1 = NTobj.edgewidth{i}(j,3);
            Y1 = NTobj.edgewidth{i}(j,4);
            X2 = NTobj.edgewidth{i}(j,5);
            Y2 = NTobj.edgewidth{i}(j,6);
            L.Position = [X1 Y1; X2 Y2];
            L.Color = 'y';
            Hlist = [Hlist L];
        end
    end
    
    display('Mark widths of desired edges. Pres Esc when done.')
        
    % Draw new interactive lines
    set(gcf,'CurrentCharacter','a');
    while true
        figure(newf)
        h = drawline('Color',[1 0.6 0.2]);
        value = double(get(gcf,'CurrentCharacter'));
        if value==27
            break;
        end
        if isempty(h)
            break;
        end          
        Hlist = [Hlist h];
    end
    
    EndAction(handles)
        
    handles.signal.String = 'Adjust width memasurements. Hit Calculate when done';
    set(handles.pushbuttonWidthCancel,'Enable','on');
    set(handles.pushbuttonWidthCalc,'Enable','on');
    % do not allow any more measurements until these are calculated or
    % cancelled
    set(handles.pushbuttonWidthMeas,'Enable','off');    
return

function pushbuttonWidthCalc_Callback(hObject, eventdata, handles)
    global NTobj Hlist
        
    handles.signal.String = 'Calculating widths...';
    for i=1:length(NTobj.edgewidth)
        NTobj.edgewidth{i} = [];
    end
    
    for i=1:length(Hlist)
        if isempty(Hlist(i))
            continue;
        end
        if norm(Hlist(i).Position(1,:) - Hlist(i).Position(2,:)) == 0
            continue;
        end
        
        [w, d, iSel] = getWidth(Hlist(i).Position);
        if (iSel == -1) % failed to find intersection
            disp('Drawn segment does not intersect the nearest edge!')
            continue;
        end
        
        if isempty(NTobj.edgewidth{iSel})
            NTobj.edgewidth{iSel} = [w d Hlist(i).Position(1,:) Hlist(i).Position(2,:)];
        else
            NTobj.edgewidth{iSel} = [NTobj.edgewidth{iSel}' [w d Hlist(i).Position(1,:) Hlist(i).Position(2,:)]']';
        end
    end
    
    for i=1:length(Hlist)
        delete(Hlist(i));
    end
    Hlist = [];
    handles.signal.String = '';
    set(handles.pushbuttonWidthCancel,'Enable','off');  
    set(handles.pushbuttonWidthCalc,'Enable','off'); 
    set(handles.pushbuttonWidthMeas,'Enable','on');
    
    disp('Width info for each measurement: width, distance along edge, 1st pt xy, 2nd pt xy')
return

function pushbuttonWidthCancel_Callback(hObject, eventdata, handles)
    global Hlist
    for i=1:length(Hlist)
        delete(Hlist(i));
    end
    Hlist = [];
    
    set(handles.pushbuttonWidthMeas,'Enable','on');
    set(handles.pushbuttonWidthCalc,'Enable','off');
    set(handles.pushbuttonWidthCancel,'Enable','off');
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Actions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function degrees = getDegrees(edgenodes, nnodes, edgeactive, selNodes)
    degrees = zeros(nnodes,1);
        E1 = edgenodes(:,1);
        E2 = edgenodes(:,2);
    for in=1:length(selNodes)
        node = selNodes(in);
        for ie=1:length(edgenodes)
            if (E1(ie)==node || E2(ie)==node) && edgeactive(ie)
                degrees(node) = degrees(node) + 1;
            end
        end
    end
return

function [nodeedges, nodenodes] = getNodeNeighbors(edgenodes, nnodes, ...
                                nodeactive, edgeactive, selNodes)
    nedges = size(edgenodes,1);
    nodeedges = zeros(nnodes,4);
    nodenodes = zeros(nnodes,4);
    ke = zeros(nnodes,1);
    kn = zeros(nnodes,1);
    for in=1:length(selNodes)
        node = selNodes(in);
        E1 = edgenodes(:,1);
        E2 = edgenodes(:,2);
        for ie = 1:nedges
            if any(edgenodes(ie,:) == node) && edgeactive(ie)
                ke(in) = ke(in) +1;
                nodeedges(node, ke(in)) = ie;
            end
            
            if E1(ie)==node && edgeactive(ie)
                kn(in) = kn(in) + 1;
                nodenodes(node, kn(in)) = E2(ie);
            end
            if E2(ie)==node && edgeactive(ie)
                kn(in) = kn(in) + 1;
                nodenodes(node, kn(in)) = E1(ie);
            end
        end
    end
return

function chain = makeNodeChains(selN, nodenodes)
% Organize all selected nodes into chains of sequential nodes    
    j = 1;
    while true
        n0 = selN(1);
        chain{j}.nodes = [n0];
        selN(selN==n0) = [];
        n = n0;
        while true
            n1 = nodenodes(n,1);
            if any(selN(:) == n1)
                chain{j}.nodes = [n1 chain{j}.nodes];
                selN(selN==n1) = [];
                n = n1;
            else
                break;
            end
        end
        
        n = n0;
        while true
            n2 = nodenodes(n,2);
            if any(selN(:) == n2)
                chain{j}.nodes = [chain{j}.nodes n2];
                selN(selN==n2) = [];
                n = n2;
            else
                break;
            end
        end
        
        if isempty(selN)
            break;
        end
        j = j+1;
    end
return

function Nodes = reorderNodes(selN, nodenodes, edgenodes)
    I = 0;
    for i=1:length(selN)
        nn = nodenodes(selN(i),1:2);
        if ~any(selN(:) == nn(1))
            Nodes = [nn(1) selN(i)];
            lastNode = nn(2);
            I = i;
            break;
        end
        if ~any(selN(:) == nn(2))
            Nodes = [nn(2) selN(i)];
            lastNode = nn(1);
            I = i;
            break;
        end
    end
    selN(I) = [];

    % Organize Nodes sequentially
    if ~isempty(selN)
        while true
            for i=1:length(selN)
                nn = nodenodes(selN(i),1:2);
                if nn(1) == Nodes(end)
                    Nodes = [Nodes selN(i)];
                    lastNode = nn(2);
                    I = i;
                    break;
                end
                if nn(2) == Nodes(end)
                    Nodes = [Nodes selN(i)];
                    lastNode = nn(1);
                    I = i;
                    break;
                end
            end
            selN(I) = [];
            
            if isempty(selN)
                break;
            end
        end
    end
    Nodes = [Nodes lastNode];
return

function [Edges, pth] = traceEdges(Nodes, edges, edgenodes, edgepath)
% Organize edges sequentially    
    nN = length(Nodes);
    nE = length(edges);
    Edges = [];
    for i=1:nN-1
        n1 = Nodes(i);
        n2 = Nodes(i+1);
        for j=1:nE
            e = edgenodes(edges(j),:);
            if e(1) == n1 && e(2) == n2 || e(1) == n2 && e(2) == n1
                Edges = [Edges edges(j)];
                if edgenodes(edges(j),:) == [n1 n2]
                    pth{i} = edgepath{edges(j)};
                else
                    pth{i} = flipud(edgepath{edges(j)});
                end
                break;
            end
        end
    end
return

function lineInd = edge2lineInd(edge, Lines, edgepath)
    lineInd = 0;
    for k=1:length(Lines)
        ex = edgepath{edge}(:,1);
        lx = Lines(k).XData(:);
        if length(lx) ~= length(ex)
            continue;
        end
        ey = edgepath{edge}(:,2);
        ly = Lines(k).YData(:);
        if nnz(ex-lx)==0 && nnz(ey-ly)==0
            lineInd = k;
            break;
        end
    end
    
    if lineInd==0 || Lines(lineInd).edgeind ~= edge
        aaa=0;
    end
    Lines(lineInd).edgeind = edge;
return

function pushbuttonMerge_Callback(hObject, eventdata, handles)
    global newf NTlocal selNodes selEdges guilock

    % merging all degree 2 nodes?
    mergeall = handles.checkboxMergeAll.Value;

    if isempty(selNodes)& ~mergeall
            return    
    end
    
    if ~NTlocal.filtered
%         % memorize selected nodes coordinates
%         memNodes = NTlocal.nodepos(selNodes,:);
%         
%         filterActiveNetwork(newf, ...
%             NTlocal.nodepos, NTlocal.edgenodes,...
%             NTlocal.nodeactive, NTlocal.edgeactive, ...
%             NTlocal.edgepath, NTobj);
%     
%         NTlocal = getNTlocal(NTobj);
%         NTlocal.filtered = true;
%         
%         % Update selNodes 
%         selN = [];
%         for in=1:size(NTlocal.nodepos,1)
%             if any(memNodes==NTlocal.nodepos(in,:))
%                 selN = [selN in];
%             end
%         end
%         selNodes = selN;
    end
    
    if (guilock)
        disp('Cannot merge, gui is locked. Finish previous operation.')
        return
    end           

    if (mergeall)
        %disp('Merging all edges is currently buggy. Do not use')
        %return
        % select all degree 2 nodes
        nnode = size(NTlocal.nodepos,1);
        degrees = getDegrees(NTlocal.edgenodes, ...
            length(NTlocal.nodeactive), NTlocal.edgeactive, 1:nnode);        
        isSel = false(1,nnode);
        isSel(selNodes) = true;
        isSel(degrees==2) = true;
        selNodes = find(isSel);
    end

    figure(newf)
    StartAction(handles, 'Merging edges...');  
    try
        degrees = getDegrees(NTlocal.edgenodes, ...
            length(NTlocal.nodeactive), NTlocal.edgeactive, selNodes);

        Lines = findobj(gca,'Type','line');
        
        % Unselect all edges
        for i=1:length(Lines)
            Lines(i).Color = 'g';
        end
        selEdges = [];

        nSel = length(selNodes);
        
        % Restrict merging only for the nodes with 2 adjoined edges,
        % unselect other nodes
        scatter = findobj(gca,'Type','scatter');
        iIgnore = [];
        for i=1:nSel
            if degrees(selNodes(i)) ~= 2
                iIgnore = [iIgnore i];
                scatter.CData(selNodes(i),:) = [1 0 0];
            end
        end
        selNodes(iIgnore) = [];
        nSel = length(selNodes);        
             
        % Edges involved
        nnodes = length(NTlocal.nodepos);
        [nodeedges, nodenodes] = getNodeNeighbors(NTlocal.edgenodes, ...
            nnodes, NTlocal.nodeactive, NTlocal.edgeactive, selNodes);
        
        edges = [];
        for i=1:nSel
            edges = [edges nodeedges(selNodes(i),1:2)];
        end
        edges = unique(edges);
        nE = length(edges);
        
        % Organize all selected nodes into chains of sequential nodes
        chain = makeNodeChains(selNodes, nodenodes);
                
        % Convert each chain into a single edge
        scatter = findobj(gca,'Type','scatter');    
        Lines = findobj(gca,'Type','line');
        for iChain=1:length(chain)
            selN = chain{iChain}.nodes;
            
            Nodes = reorderNodes(selN, nodenodes);            
            [Edges pth] = ...
                traceEdges(Nodes, edges, NTlocal.edgenodes,...
                NTlocal.edgepath);
            NTlocal.edgenodes(Edges(1),:) = ...
                [Nodes(1) Nodes(end)];

            newpath = pth{1};
            for i=2:length(pth)
                newpath = [newpath(1:end-1,:); pth{i}];
            end            
            
            for i=1:length(selN)
                scatter.SizeData(selN(i)) = 0.001;
                scatter.CData(selN(i),:) = [0 0 0];
                edge1 = nodeedges(selN(i),1);
                edge2 = nodeedges(selN(i),2);
                
                lineInd1 = edge2lineInd(edge1, Lines, NTlocal.edgepath);
                lineInd2 = edge2lineInd(edge2, Lines, NTlocal.edgepath);
                if i>1
                    Lines(lineInd1).Visible = 'off';
                end
                if i==1
                    Lines(lineInd1).XData = newpath(:,1);
                    try
                        Lines(lineInd1).YData = newpath(:,2);
                    catch
                        aaa=0;
                    end
                end
                Lines(lineInd2).Visible = 'off';
                
                if i>1
                    NTlocal.edgeactive(edge1) = false;
                end
                NTlocal.edgeactive(edge2) = false;
            end
            NTlocal.nodeactive(selN) = false;
    
            NTlocal.edgepath{Edges(1)} = newpath;            
        end
        
        % clear selected nodes 
        selNodes = [];
        % local network topology has changed
        NTlocal.filtered = 0;

    catch exception
        guilock = false;
        disp(getReport(exception))
    end
    EndAction(handles)
return

function removeSelected()
    global NTlocal newf selNodes selEdges 

    if isempty(selNodes) & isempty(selEdges)
        return
    end
        
    figure(newf);
    % Inactivate selected nodes
    scatter = findobj(gca,'Type','scatter');    
    scatter.SizeData(selNodes) = 0.001;
    NTlocal.nodeactive(selNodes) = false; 

    % Add adjacent edges to selected for later deletion
    for j=1:length(selNodes)
        k = selNodes(j);
        scatter.CData(k,:) = [0 0 0];
        for ec=1:length(NTlocal.edgepath)
            if any(NTlocal.edgenodes(ec,:) == k)
                selEdges = [selEdges ec]; 
            end
        end
    end
    selEdges = sort(unique(selEdges));
    
    % Inactivate selected edges
    Lines = findobj(gca,'Type','line');
    ind = [];
    for j=1:length(Lines)
        if any(selEdges(:) == Lines(j).edgeind) 
            Lines(j).Visible = 'off';
        end
    end
    NTlocal.edgeactive(selEdges) = false; 
    
    selNodes = [];
    selEdges = [];
    NTlocal.filtered = false;
return

function pushbuttonRemoveSelected_Callback(hObject, eventdata, handles)
    StartAction(handles, 'Removing selected...')
    removeSelected();
    EndAction(handles)
return

function pushbuttonFilterActive_Callback(hObject, eventdata, handles)
    global newf NTobj NTlocal edgeplotH guilock
    
    % convert the current active network in the GUI into the workspace
    % object

    %if ~(NTlocal.filtered)

     filterActiveNetwork(newf, NTlocal.nodepos, NTlocal.edgenodes,...
            NTlocal.nodeactive, NTlocal.edgeactive, NTlocal.edgepath,...
            NTlocal.edgewidth,NTobj);

     % also replot the network in the object
     NTlocal = getNTlocal(NTobj);
    %end
    dispNetWithImage() 
return

function mapall2act = filterActiveNetwork(newf, nodepos, edgenodes,...
                                  nodeact, edgeact, edgepath, edgewidth, NT)
    % from a full list of nodes and edges (some inactive)
    % filter out only the active ones and make a clean network obj

    %NT = NetworkObj();

    prevrootnode = NT.rootnode;

    % copy over active node positions
    NT.nodepos = nodepos(nodeact,:);
    
    
    % copy over active edges    
    newedgenodes = edgenodes(edgeact,:);
    NT.edgenodes = zeros(nnz(edgeact),2);
    
    % Clean up which nodes the edges are pointing to:
    
    % this maps the list of all nodes to their index in the list of active nodes
    mapall2act = zeros(1,size(nodepos,1));
    mapall2act(nodeact) = 1:nnz(nodeact);
    
    NT.edgenodes(:,1) = mapall2act(newedgenodes(:,1));
    NT.edgenodes(:,2) = mapall2act(newedgenodes(:,2));
    
    % copy over active edge paths
    NT.edgepath = edgepath(edgeact);

    % copy over edge widths  
    NT.edgewidth = edgewidth(edgeact);
    
    % set up all the other arrays in the network
    NT.setupNetwork();
    
    % deal with cumulative and total edge lengths
    NT.setCumEdgeLen(1:NT.nedge, true);

    % update root node
    if (~isnan(prevrootnode))
        NT.rootnode = mapall2act(prevrootnode);

        % update directionality info
        isset = false(NT.nedge,1);
        wasreversed = false(NT.nedge,1);
        [isset,wasreversed] = directedTreeEdges(NT,NT.rootnode,isset,wasreversed);
    end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   UTILS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PickableOn(obj)
    global guilock newf   
    guilock = true;
    figure(newf)
    for i=1:length(obj)
        if(~isempty(obj(i)))
            obj(i).PickableParts = 'all';
            obj(i).HitTest = 'on';
        end
    end
return

function PickableOff(obj)
    global guilock newf   
    guilock = false;
    figure(newf)
    set(gcf,'Pointer','arrow');
    for i=1:length(obj)
        if(~isempty(obj(i)))
            obj(i).PickableParts = 'none';
            obj(i).HitTest = 'off';
        end
    end
return

function ActionsEnable(handles)
    hb = findobj('Style','pushbutton');
    for i=1:length(hb)
        hb(i).Enable = 'On';
    end
    handles.sliderContrast.Enable = 'On';
    handles.sliderBrightness.Enable = 'On';
return

function ActionsDisable(handles)
    hb = findobj('Style','pushbutton');
    for i=1:length(hb)
        hb(i).Enable = 'Off';
    end
    handles.sliderContrast.Enable = 'Off';
    handles.sliderBrightness.Enable = 'Off';
return

function StartAction(handles, str)
global newf guilock
    ActionsDisable(handles);
    guilock = true;
    handles.signal.String = str;
    figure(newf)
%     set(gcf,'Pointer','watch');
return

function EndAction(handles)
global newf guilock
    ActionsEnable(handles);
    guilock = false;
    handles.signal.String = '';
    figure(newf)
    set(gcf,'Pointer','Arrow');
return


% --- Executes on button press in pushbuttonColorByDegree.
function pushbuttonColorByDegree_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonColorByDegree (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global NTlocal newf NTobj

nedge = size(NTlocal.edgenodes,1);
nnode = size(NTlocal.nodepos,1);
degs = zeros(nnode,1);
for ec = 1:nedge
    if (NTlocal.edgeactive(ec))
        n1 = NTlocal.edgenodes(ec,1);
        n2 = NTlocal.edgenodes(ec,2);
        if (NTlocal.nodeactive(n1)); degs(n1) = degs(n1)+1; end
        if (NTlocal.nodeactive(n2)); degs(n2) = degs(n2)+1; end
    end
end

colors = lines(max(degs));
colors(1,:) = colors(1,:)*1.5;

figure(newf)
scatter = findobj(gca,'Type','scatter');

for i=1:nnode
    if (degs(i)>0)
        scatter.CData(i,:) = colors(degs(i),:);
    end
end

return

% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over pushbuttonColorByDegree.
function pushbuttonColorByDegree_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to pushbuttonColorByDegree (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in ShowLoops.
function ShowLoops_Callback(hObject, eventdata, handles)
% hObject    handle to ShowLoops (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global newf NTlocal NTobj

% make full network object
NTcur = copy(NTobj);
filterActiveNetwork(newf, NTlocal.nodepos, NTlocal.edgenodes,...
    NTlocal.nodeactive, NTlocal.edgeactive, NTlocal.edgepath,...
    NTcur);

% set up a connectivity matrix
nnode = size(NTlocal.nodepos,1);
nedge = size(NTlocal.edgenodes,1);
A = zeros(nnode,nnode);
for ec = 1:nedge
    if (NTlocal.edgeactive(ec,:))
        n1 = NTlocal.edgenodes(ec,1); n2 = NTlocal.edgenodes(ec,2);
        A(n1,n2) = 1;
        A(n2,n1) = 1;
    end
end
% make into a graph structure
NTgraph = graph(A);

cycles = cyclebasis(NTgraph);
cycleNodes = horzcat(cycles{1:end});

figure(newf)
scatter = findobj(gca,'Type','scatter');

for i=cycleNodes
    scatter.CData(i,:) = [0,0,1];
end

return


% --- Executes on button press in checkboxMergeAll.
function checkboxMergeAll_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxMergeAll (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxMergeAll


% --- Executes on button press in pushbuttonSetParentNode.
function pushbuttonSetParentNode_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonSetParentNode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global NTobj NTlocal selNodes newf


% set the parent node to the (single) selected node
if (isempty(selNodes))
    handles.signal.String = 'No nodes selected. Cannot set parent node.';
    return
elseif (length(selNodes)>1)
    handles.signal.String = 'Multiple nodes selected. Cannot set parent node.';
    return
end


% Update the original network structure. 
locrootnode = selNodes(1); % local index of root node

mapall2act = filterActiveNetwork(newf, NTlocal.nodepos, NTlocal.edgenodes,...
        NTlocal.nodeactive, NTlocal.edgeactive, NTlocal.edgepath,...
        NTobj);

% new index of root node
rootnode = mapall2act(locrootnode);
handles.signal.String = sprintf('Set root node in original network object: %d',rootnode);
NTobj.rootnode = rootnode;
selNodes = [];

NTlocal = getNTlocal(NTobj);   
dispNetWithImage() 

% reverse edge directions to make network into a directed tree
% starting from the root node
isset = false(NTobj.nedge,1);
wasreversed = false(NTobj.nedge,1);
[isset,wasreversed] = directedTreeEdges(NTobj,rootnode,isset,wasreversed);
handles.signal.String = sprintf('Set root node in original network object: %d. Set directed edges.',rootnode);
return


% --- Executes on button press in pushbuttonReplotNetwork.
function pushbuttonReplotNetwork_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonReplotNetwork (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global NTobj NTlocal 

% replot the network in the GUI based on the current workspace object (in
% NTobj)

% also replot the network in the object
NTlocal = getNTlocal(NTobj);
%end
dispNetWithImage()
