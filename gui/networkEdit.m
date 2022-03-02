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

% Last Modified by GUIDE v2.5 30-Jan-2022 10:52:06

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
    % index mapping from the current network to the plotted elements
    global nodenet2plot edgenet2plot
    
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
    if isempty(NTobj.edgewidth)
        NTobj.edgewidth = cell(NT.nedge,1);
    end
    
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
global newf NTobj imgObj plotoptObj nodeplotH edgeplotH imageH selNodes nodenet2plot edgenet2plot
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
    
    % initialize mapping of network nodes/edges to plotted elements
    nodenet2plot = 1:NTobj.nnode;
    edgenet2plot = 1:NTobj.nedge;
return

function plotNet()
    global newf NTobj plotoptObj nodeplotH edgeplotH nodenet2plot edgenet2plot
    
    newf=figure(1);
    
    scatter = findobj(gca,'Type','scatter');
    delete(scatter);
    
    Lines = findobj(gca,'Type','line');
    delete(Lines);

    hold on   
    plotoptObj.datatipindex = true;
    [nodeplotH,edgeplotH] = NTobj.plotNetwork(plotoptObj);
    hold off
    
    nodenet2plot = 1:NTobj.nnode;
    edgenet2plot = 1:NTobj.nedge;
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
    global guilock NTobj newf selNodes
    
    if (guilock)
        disp('Cannot select nodes, gui is locked. Finish previous operation.')
        return
    end
    StartAction(handles, "Pick points to select desired nodes. Hit Esc when done.")
    %set(gcf,'Pointer','Arrow');
    %ind = selectNode(true, [0 0 1]);   
    
    figure(newf)
    
    [P,H]= selectMultiplePoints();
    
    % get all the roi points, including the ones drawn previously
    H = findobj(gca,'Type','images.roi.Point');
    goodpt = false(length(H),1);
    for pc = 1:length(H)
        goodpt(pc) = ~isempty(H(pc).Position);
    end
    H = H(goodpt);    
    P = vertcat(H.Position);
    
    % find the nodes closest to the selected points
    if (~isempty(P))
    selNodes = knnsearch(NTobj.nodepos,P);
    
    % snap draggable node positions to the network
    for pc = 1:length(H)
        H(pc).Position = NTobj.nodepos(selNodes(pc),:);
    end
    end
    
    
    EndAction(handles)
return

function pushbuttonUnselectAllNodes_Callback(hObject, eventdata, handles)
    global newf selNodes
        
    figure(newf)
    nodeH = findobj(gca,'Type','images.roi.Point');
    delete(nodeH)
    
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
    
    % make roi points on the selected nodes
    selNodes = find(in);
    
    for pc = 1:length(selNodes)
        drawpoint('Position',NTobj.nodepos(selNodes(pc),:));
    end
    
    hold off
return

function pushbuttonAddNode_Callback(hObject, eventdata, handles)
    global newf NTobj nodeplotH edgeplotH nodenet2plot
    
    StartAction(handles, ...
        'Click on positions of desired nodes. To finish hit Esc.')
    try
        [P,H]= selectMultiplePoints();
        
        set(gcf,'Pointer','arrow');
       
        for i = 1:length(P)
            connect{i} = [];
            delete(H(i));
        end
        NTobj.addNodes(P,connect)
                
        handles.signal.String = 'WAIT...';
        %plotNet();
        
        % replot scatter points only
        delete(nodeplotH);
        hold all
        cols = [ones(NTobj.nnode,1) zeros(NTobj.nnode,2)];
        nodeplotH = scatter(NTobj.nodepos(:,1),NTobj.nodepos(:,2),20*ones(1,NTobj.nnode),cols,'filled');
        hold off
        nodenet2plot = 1:NTobj.nnode;
        
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
        disp(getReport(exception))
    end
    EndAction(handles)
return

function [P,H]= selectMultiplePoints()
    % select multiple points
    % for node selection or adding new nodes
    % returns position coords in P and handles in H    
    try
        P = []; H=[];
        
        i=0;
        while true
            h = drawpoint;
            
            if isempty(h.Position)
                break;
            end  
            i = i+1;
            P(i,:) = h.Position;
            H = [H h];
        end
        set(gcf,'Pointer','arrow');
    catch exception
        disp(getReport(exception))
    end    
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Manage edges
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function iSel = selectEdge(color)
    global newf NTobj selEdges edgeplotH guilock selectFirst
    
    try
        figure(newf)
        iSel=[];  
      
        display('Select points marking desired edges. Press escape when done.')
        
        [P,H] = selectMultiplePoints();
        
        % figure which edges are nearest selected points
        for pc = 1:size(P,1)
            [mindist,minec,minfrac,minpt] = NTobj.getNearestEdge(P(pc,:));
            
            iSel(end+1) = minec;
                        
            % remove the original drawn point
            delete(H(pc))
        end
        
    catch exception
        disp(getReport(exception))
    end
    
return

function pushbuttonSelectEdge_Callback(hObject, eventdata, handles)
    global guilock selEdges NTobj
    if (guilock)
        disp('Cannot select edges, gui is locked. Finish previous operation.')
        return
    end
    StartAction(handles, "Click on (or near) desired edges to select. Hit Esc to finish.")
    set(gcf,'Pointer','Arrow');
    
    isselected = false(1,NTobj.nedge);
    isselected(selEdges) = true;
    
    iSel = selectEdge('b');
    % flip selection state of the picked edges (either selects or
    % unselects)     
    isselected(iSel) = ~isselected(iSel);
    selEdges = find(isselected);
    
    % get rid of prior polylines
    tmp= findobj(gca,'Type','images.roi.PolyLine');
    delete(tmp);
    
    % mark selected edges with roi polyline
    for ec = selEdges
        drawpolyline('Position',NTobj.edgepath{ec},'Color','b','InteractionsAllowed','none');
    end
    
    EndAction(handles)
return


function pushbuttonUnselectAllEdges_Callback(hObject, eventdata, handles)
    global newf selEdges

    figure(newf)
    L = findobj(gca,'Type','images.roi.PolyLine');
    delete(L);
    
    selEdges = [];
return

function pushbuttonAddEdge_Callback(hObject, eventdata, handles)
    global newf NTobj edgeplotH guilock edgenet2plot
    
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
        
        edgenet2plot(NTobj.nedge) = max(edgenet2plot)+1;

       % CLF = hypot(diff(x), diff(y));   % Calculate integrand from x,y derivatives
       % NTobj.edgelens(i) = trapz(CLF); % Integrate to calculate arc length

        NTobj.edgeedges(i,:,:) = NTobj.edgeedges(i-1,:,:);

        hold on
        edgeplotH(i) = plot(x,y, 'Color','g', 'LineWidth',2);
        hold off
        edgeplotH(i).addprop('edgeind');
        edgeplotH(i).edgeind = i;
        
%         dttemplate = edgeplotH(i).DataTipTemplate;
%         dttemplate.FontSize=6;
%         dttemplate.DataTipRows(1).Value = i*ones(size(NTobj.edgepath{i},1),1);
%         dttemplate.DataTipRows(1).Label = '';
%         dttemplate.DataTipRows(2:end) = [];
%         dttemplateset = true;
        
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

function iSel = findNearestEdge(xy)
    global NTobj

    dMin = inf;
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
% merge at selected nodes. This will trigger redrawing of full network
    global NTobj newf selNodes selEdges nodeplotH guilock edgenet2plot nodenet2plot

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
        L = findobj(gca,'Type','images.roi.PolyLine');
        delete(L);    
        selEdges = [];
        
        
        % get selected node handles
        nodeH = findobj(gca, 'Type','images.roi.Point');       
    
        if (isempty(nodeH))
            disp('No nodes selected for merging.')
            return
        end
    
        % positions of points
        selnodepos = vertcat(nodeH.Position);
    
        if (isempty(selnodepos))
            disp('No nodes selected for merging.')
            delete(nodeH)
            return
        else
            % get selected network nodes
            selNodes = knnsearch(NTobj.nodepos,selnodepos);
        end
        
        deg2ind = find(NTobj.degrees(selNodes)==2);
        selNodes = selNodes(deg2ind);
        
        if (isempty(selNodes))
            disp('No deg 2 nodes selected')
            return
        end
        
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
        
        % delete all picked nodes
        delete(nodeH)        
        
        figure(newf)
        % replot network
        plotNet();
    catch exception
        disp(getReport(exception))
    end
    EndAction(handles)
return

function removeSelected()
    global NTobj newf selEdges nodeplotH edgeplotH guilock nodenet2plot edgenet2plot

    figure(newf)    
    
    %% check node numbering in plot and network
    for pc = 1:NTobj.nnode
        ind = nodenet2plot(pc);
        diffval = norm(NTobj.nodepos(pc,:) - [nodeplotH.XData(ind), nodeplotH.YData(ind)]);
        if (diffval>1e-10)
            disp([pc diffval])
            error('Network node data and plotted scatter plot indexing have become mismatched. Something is broken. ')
        end
    end
    
    
    % selected node handles
    nodeH = findobj(gca, 'Type','images.roi.Point');       
    
    if (isempty(nodeH) & isempty(selEdges))
        disp('No nodes or edges selected to remove')
        return
    end
    
    if (isempty(nodeH))
        selnodepos = [];
    else
        selnodepos = vertcat(nodeH.Position);
    end
    
%     figure(newf);
%     guilock = true;
%     set(gcf,'Pointer','watch');
%     handles.signal.String = 'Removing selected...';
   % try
        % remove nodes and adjacent edges to those nodes
        if (size(selnodepos,1)>0)         
             selNodes = knnsearch(NTobj.nodepos,selnodepos);
    
            nedge0 = NTobj.nedge;
            stillthere = false(1,nedge0);
            
            
            dokeep = true(1,NTobj.nnode);
            dokeep(selNodes) = false;
            keepind = find(dokeep);
            mapold2newedge = zeros(1,NTobj.nedge);
            [mapold2new,mapnew2oldedge] = NTobj.keepNodes(keepind);

            tmp = mapold2new(keepind);            
            mapnew2oldnode(tmp) = keepind;
            
            % make removed nodes invisible
            nodeplotH.SizeData(nodenet2plot(selNodes)) = nan;
            
            % which edges are still there, and which got removed
            stillthere(mapnew2oldedge) = true;
            rmedges = find(~stillthere);
            
            % make removed edges not show
            for ec = 1:length(rmedges)
                edgeplotH(edgenet2plot(rmedges(ec))).LineStyle = 'none';
            end
            
            % update mapping from new network node/edges to plot elements
            edgenet2plot = edgenet2plot(mapnew2oldedge);
            nodenet2plot = nodenet2plot(mapnew2oldnode);
            
            % update index of selected edges
            mapold2newedge(mapnew2oldedge) = 1:NTobj.nedge;
            selEdges = mapold2newedge(selEdges);
            selEdges = selEdges(selEdges>0);
        end

        % remove additional edges 
        dokeep = true(1,NTobj.nedge);
        dokeep(selEdges) = false;
        keepind = find(dokeep);
        mapnew2oldedge = NTobj.keepEdges(keepind);
       
        % make removed edges not show
        for ec = 1:length(selEdges)
            edgeplotH(edgenet2plot(selEdges(ec))).LineStyle = 'none';
        end
                
        selNodes = [];
        selEdges = [];
        
        edgenet2plot = edgenet2plot(mapnew2oldedge);
        
   % catch exception
   %     disp(getReport(exception))
   % end
    
    % delete marked node handles
    delete(nodeH);
    
    % delete selected edge handles
    edgeH = findobj(gca,'Type','images.roi.PolyLine');
    delete(edgeH);
    
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
