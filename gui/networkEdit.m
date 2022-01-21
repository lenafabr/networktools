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

% Last Modified by GUIDE v2.5 27-Sep-2021 11:39:00

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
    global newf NTobj imgObj plotoptObj selNodes selEdges guilock
    global selectFirst
% Choose default command line output for networkEdit
handles.output = hObject;
selectFirst = true;


% Update handles structure
guidata(hObject, handles);

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

    newf =[];
    selNodes = [];
    selEdges = [];
    if (~isempty(NTobj))
        NTobj.edgewidth = cell(NTobj.nedge,1);
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

    % UIWAIT makes networkEdit wait for user response (see UIRESUME)
% uiwait(handles.mainFig);

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
    global newf NTobj imgObj plotoptObj selNodes selEdges fileName
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
    NTobj.edgewidth = cell(NT.nedge,1);
        
    figure(1);       
    imageH = imshow(imgObj,[]);
%     imgCData0 = imageH.CData;

    dispNetWithImage();
    figure(newf);
return

function menuLoad_Callback(hObject, eventdata, handles)
    LoadData();
    ActionsEnable(handles);
return

function menuSave_Callback(hObject, eventdata, handles)
    global NTobj imgObj plotoptObj fileName
    
    NT = NTobj;
    img = imgObj;
    plotopt = plotoptObj;
    newFile = [fileName(1:end-4) '_modified.mat'];
    uisave({'NT','img','plotopt'}, newFile);
return

function menuClear_Callback(hObject, eventdata, handles)
    global newf NTobj selNodes 
    
    close(figure(1))
    newf =[];
    NTobj = [];
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
            imageexists = isvalid(imageH);
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

function plotNet()
    global newf NTobj plotoptObj nodeplotH edgeplotH
    
    newf=figure(1);
    
    scatter = findobj(gca,'Type','scatter');
    delete(scatter);
    
    Lines = findobj(gca,'Type','line');
    delete(Lines);

    hold on   
    plotoptObj.datatipindex = true;
    [nodeplotH,edgeplotH] = NTobj.plotNetwork(plotoptObj);
    hold off
return

function sliderContrast_Callback(hObject, eventdata, handles)
    global imageH  imgObj
    try
        a = get(hObject,'Value');
        hSlider2 = findobj('Tag', 'sliderBrightness');
        b = hSlider2.Value;
        imageH.CData = imgObj * a + b;
              
        imageH.CData(imageH.CData>1) = 1;
        imageH.CData(imageH.CData<0) = 0;
    catch exception
        disp(getReport(exception))
    end
return

function sliderBrightness_Callback(hObject, eventdata, handles)
    global imageH B imgObj
    try
        hSlider1 = findobj('Tag', 'sliderContrast');
        a = hSlider1.Value;
        b = get(hObject,'Value');
        imageH.CData = imgObj*a+b
        
        imageH.CData(imageH.CData>1) = 1;
        imageH.CData(imageH.CData<0) = 0;
    catch exception
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Manage nodes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ind = selectNode(addSelected, color)
    global newf nodeplotH selNodes guilock selectFirst

    disp('Use datatips to select desired nodes. Then hit any keyboard key (while the figure window is active).')
   
    try
        figure(newf)
%         if selectFirst
%             imitateEnter()
%         end

        scatter = findobj(gca,'Type','scatter');
        PickableOn(nodeplotH, scatter);
        
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
            PickableOff(nodeplotH, scatter);
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
        disp(getReport(exception))
    end
    if ~exist('newf','var')
        return;
    end
    figure(newf)
    scatter = findobj(gca,'Type','scatter');
    PickableOff(nodeplotH, scatter);
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
    global newf NTobj selNodes selEdges
    
    selEdges = unique(selEdges);
    figure(newf);
    hold on
    
    roi = drawpolyline;
    xp = roi.Position(:,1);
    yp = roi.Position(:,2);
    xc = [xp' roi.Position(1,1)]';
    yc = [yp' roi.Position(1,2)]';
    %plot(xp,yp, 'm')
    
    x = NTobj.nodepos(:,1);
    y = NTobj.nodepos(:,2);
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
    global newf NTobj selNodes
    
    try
        figure(newf);

        roi = drawpoint;
        xy = roi.Position;
        delete(roi);
        set(gcf,'Pointer','watch');
        handles.signal.String = 'WAIT...';

        NTobj.nnode = NTobj.nnode + 1;
        nnode = NTobj.nnode;
        NTobj.nodepos(nnode,:) = xy';
        NTobj.degrees(nnode) = 0;
        NTobj.nodenodes(nnode,:) = [0 0 0 0];
        NTobj.nodeedges(nnode,:) = [0 0 0 0];

        plotNet();
    catch exception
        disp(getReport(exception))
    end
    figure(newf)
    set(gcf,'Pointer','arrow');
    handles.signal.String = '';
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Manage edges
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function iSel = selectEdge(color)
    global newf NTobj selEdges edgeplotH guilock selectFirst
    
    L = findobj(gca,'Type','line');
    PickableOn(edgeplotH, L);
    try
        figure(newf)
%         if selectFirst
%             imitateEnter()
%         end

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
            PickableOff(edgeplotH, L);
            return
        end

        if ~isempty(datatips)
            for idt=1:length(datatips)
                Ln = datatips(idt).Parent;
                Ln.Color = color;
                iSel = [iSel Ln.edgeind];
            end

            delete(datatips)
            selEdges = [selEdges iSel];
        end
    catch exception
        disp(getReport(exception))
    end
    figure(newf)
    L = findobj(gca,'Type','line');
    PickableOff(edgeplotH, L);
return

function pushbuttonSelectEdge_Callback(hObject, eventdata, handles)
    global guilock
    if (guilock)
        disp('Cannot select edges, gui is locked. Finish previous operation.')
        return
    end
    StartAction(handles, "Use datatips to select edges. Then hit Enter to complete, Esc to cancel")
    set(gcf,'Pointer','Arrow');
    iSel = selectEdge('b');
    EndAction(handles)
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
    global newf NTobj edgeplotH guilock
    
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

        NTobj.setCumEdgeLen(i,true);
        NTobj.edgevals{i} = [];
        NTobj.edgewidth{i} = [];

       % CLF = hypot(diff(x), diff(y));   % Calculate integrand from x,y derivatives
       % NTobj.edgelens(i) = trapz(CLF); % Integrate to calculate arc length

        NTobj.edgeedges(i,:,:) = NTobj.edgeedges(i-1,:,:);

        hold on
        edgeplotH(i) = plot(x,y, 'Color','g', 'LineWidth',2);
        hold off
        edgeplotH(i).addprop('edgeind');
        edgeplotH(i).edgeind = i;
    catch exception
        disp(getReport(exception))
    end
    guilock = false;
return
 
function ind = findNearestNode(nodepos, xy)
    dNodes = nodepos - xy;
    d2 = dNodes(:,1).^2 + dNodes(:,2).^2;
    M = min(d2);
    ind = find(d2 == M);       
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Actions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function chain = makeNodeChains(selN)
% Organize all selected nodes into chains of sequential nodes
    global NTobj selNodes
    
    j = 1;
    while true
        n0 = selN(1);
        chain{j}.nodes = [n0];
        selN(selN==n0) = [];
        n = n0;
        while true
            n1 = NTobj.nodenodes(n,1);
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
            n2 = NTobj.nodenodes(n,2);
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
    return;
return

function Nodes = reorderNodes(selN)
    global NTobj
    
    I = 0;
    for i=1:length(selN)
        nn = NTobj.nodenodes(selN(i),1:2);
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
                nn = NTobj.nodenodes(selN(i),1:2);
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

function [Edges pth] = traceEdges(Nodes, edges)
% Organize edges sequentially
    global NTobj
    
    nN = length(Nodes);
    nE = length(edges);
    Edges = [];
    for i=1:nN-1
        n1 = Nodes(i);
        n2 = Nodes(i+1);
        for j=1:nE
            e = NTobj.edgenodes(edges(j),:);
            if e(1) == n1 && e(2) == n2 || e(1) == n2 && e(2) == n1
                Edges = [Edges edges(j)];
                if NTobj.edgenodes(edges(j),:) == [n1 n2]
                    pth{i} = NTobj.edgepath{edges(j)};
                else
                    pth{i} = flipud(NTobj.edgepath{edges(j)});
                end
                break;
            end
        end
    end
return

function pushbuttonMerge_Callback(hObject, eventdata, handles)
    global NTobj newf selNodes selEdges nodeplotH guilock

    if isempty(selNodes)
        return
    end

    if (guilock)
        disp('Cannot merge, gui is locked. Finish previous operation.')
        return
    end
    
    figure(newf)
    StartAction(handles, 'Merging edges...');  
    try
        % Unselect all edges
        L = findobj(gca,'Type','line');
        for i=1:length(L)
            L(i).Color = 'g';
        end
        selEdges = [];

        nSel = length(selNodes);
        
        % Restrict merging only for the nodes with 2 adjoined edges,
        % unselect other nodes
        scatter = findobj(gca,'Type','scatter');
        iIgnore = [];
        for i=1:nSel
            if NTobj.degrees(selNodes(i)) ~= 2
                iIgnore = [iIgnore i];
                scatter.CData(selNodes(i),:) = [1 0 0];
            end
        end
        selNodes(iIgnore) = [];
        nSel = length(selNodes);

        % Edges involved
        edges = [];
        for i=1:nSel
            edges = [edges NTobj.nodeedges(selNodes(i),1:2)];
        end
        edges = unique(edges);
        nE = length(edges);

        % Organize all selected nodes into chains of sequential nodes
        chain = makeNodeChains(selNodes);
        
        dokeepEdge = true(1,NTobj.nedge);
        dokeepNode = true(1,NTobj.nnode);
        
        % Convert each chain into a single edge
        for iChain=1:length(chain)
            selN = chain{iChain}.nodes;
            dokeepNode(selN) = false;
            
            Nodes = reorderNodes(selN);
            nN = length(Nodes);
            
            [Edges pth] = traceEdges(Nodes, edges);
            NTobj.edgenodes(Edges(1),:) = [Nodes(1) Nodes(end)];

            newpath = pth{1};
            for i=2:length(pth)
                newpath = [newpath(1:end-1,:); pth{i}];
            end
            NTobj.edgepath{Edges(1)} = newpath;
            
            dokeepEdge(Edges(2:end)) = false;
        end
        
        keepindEdge = find(dokeepEdge);
        keepindNode = find(dokeepNode);

        NTobj.setupNetwork()
        NTobj.keepEdges(keepindEdge);
        NTobj.keepNodes(keepindNode);
        
        figure(newf)
        plotNet();
    catch exception
        disp(getReport(exception))
    end
    EndAction(handles)
return

function removeSelected()
    global NTobj newf selNodes selEdges nodeplotH edgeplotH guilock

    if isempty(selNodes) & isempty(selEdges)
        return
    end
    
%     figure(newf);
%     guilock = true;
%     set(gcf,'Pointer','watch');
%     handles.signal.String = 'Removing selected...';
    try
        % remove nodes and adjacent edges to those nodes
        if (~isempty(selNodes))
            dokeep = true(1,NTobj.nnode);
            dokeep(selNodes) = false;
            keepind = find(dokeep);
            mapold2newedge = zeros(1,NTobj.nedge);
            [~,mapnew2oldedge] = NTobj.keepNodes(keepind);

            % update index of selected edges
            mapold2newedge(mapnew2oldedge) = 1:NTobj.nedge;
            selEdges = mapold2newedge(selEdges);
        end

        % remove additional edges 
        dokeep = true(1,NTobj.nedge);
        dokeep(selEdges) = false;
        keepind = find(dokeep);
        NTobj.keepEdges(keepind);

        selNodes = [];
        selEdges = [];
        plotNet(); 
    catch exception
        disp(getReport(exception))
    end
    
%     handles.signal.String = '';
%     guilock = false;
%     figure(newf)
%     set(gcf,'Pointer','arrow');
return

function pushbuttonRemoveSelected_Callback(hObject, eventdata, handles)
    StartAction(handles, 'Removing selected...')
    removeSelected();
    EndAction(handles)
return

function iSel = findNearestEdge(xy)
    global NTobj

    dMin = 1E20;
    iSel =-1;
    for i=1:NTobj.nedge
        d = NTobj.edgepath{i} - xy;
        d2 = d(:,1).^2 + d(:,2).^2;
        m = min(d2);
        if m<dMin
            iSel = i;
            dMin = m;
        end
    end
return

function pushbuttonEdgeWidths_Callback(hObject, eventdata, handles)
    global newf NTobj guilock
              
    %warning('Width measurement has not been tested, is probably buggy, and should not be used.')
    %return
    
    if (guilock)
        disp('Cannot merge, gui is locked. Finish previous operation.')
        ind = [];
        return
    end
    
%     set(gcf,'Pointer','watch');
%     handles.signal.String = 'WAIT...';
    StartAction(handles, 'Calculating edge width...')

    figure(newf);
    hold on
    try
        h = drawline();
        ep = h.Position;
        X1 = ep(1,1); Y1 = ep(1,2); X2 = ep(2,1); Y2 = ep(2,2);
        xy = [0.5*(X1+X2) 0.5*(Y1+Y2)];
        iSel = findNearestEdge(xy);

        edge = NTobj.edgepath{iSel};   
        for j=1:size(edge,1)-1
            x1 = edge(j,1);
            y1 = edge(j,2);
            x2 = edge(j+1,1);
            y2 = edge(j+1,2);

            p1 = polyfit([x1 x2], [y1 y2], 1);
            p2 = polyfit([X1 X2], [Y1 Y2],1);

            if abs(x2-x1)<1E-10
                xCros = x1;
                yCros = Y1 + (xCros - X1)*(Y2 - Y1)/(X2 - X1);
            else
                if abs(y1-y2)<1E-10
                    yCros = y1;
                    xCros = X1 + (yCros - Y1)*(X2 - X1)/(Y2 - Y1);
                else
                    xCros = fzero(@(x) polyval(p1-p2,x),3);
                    yCros = polyval(p1,xCros);
                end
            end
            if (xCros-x1)*(xCros-x2)<=0 & (yCros-y1)*(yCros-y2)<=0
                break;
            end

            plot(xCros, yCros, 'r+', 'LineWidth', 2)
        end
        hold off

        w = norm([X2-X1, Y2-Y1]);
        p = NTobj.edgepath{iSel};
        cum = NTobj.cumedgelen{iSel};
        d = cum(j) + norm([xCros yCros] - p(j,:));
        if isempty(NTobj.edgewidth{iSel})
            NTobj.edgewidth{iSel} = [w d];
        else
            NTobj.edgewidth{iSel} = [NTobj.edgewidth{iSel}' [w d]']';
        end
%         figure(newf);
    catch exception
        disp(getReport(exception))
%         figure(newf);
%         set(gcf,'Pointer','arrow');
    end
    
%     set(gcf,'Pointer','arrow');
%     handles.signal.String = '';
    EndAction(handles)
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   UTILS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PickableOn(obj1, obj2, on)
    n = length(obj1);
    if length(obj1) ~= n
        return;
    end
    
    guilock = true;
    for i=1:n
        obj1(i).PickableParts = 'all';
        obj2(i).PickableParts = 'all';
        obj1(i).HitTest = 'on';
        obj2(i).HitTest = 'on';
    end
return

function PickableOff(obj1, obj2, on)
    global guilock newf
        
    n = length(obj1);
    if length(obj1) ~= n
        return;
    end
    
    guilock = false;
    figure(newf)
    set(gcf,'Pointer','arrow');
    for i=1:n
        obj1(i).PickableParts = 'none';
        obj2(i).PickableParts = 'none';
        obj1(i).HitTest = 'off';
        obj2(i).HitTest = 'off';
    end
return

function imitateEnter()
% Imitation of clicking Enter on the figure
    global selectFirst
    try
        P = get(gcf, 'Position');
        %             SS = get(0,'screensize');
        hroot=groot;
        set(hroot,'PointerLocation',[P(1)+10, P(2)+10]);
        pause(2);
        inputemu('key_normal','\ENTER');
    catch exception
        disp(getReport(exception))
        selectFirst = false;
        return
    end
    selectFirst = false;
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
