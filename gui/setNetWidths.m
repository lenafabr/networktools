function [NT,saveedgewidth] = setNetWidths(NT,imgfilename)
%% function to set network edgewidth (diameters) through gui
% Later: add the option of deleting nodes and corresponding edges
%% Function has 2 options: 
% Manual entry (click on edge, set diameter from other measurements)
% Measure tool (measure directly on matlab)
% Measure tool stores edgewidth as NT.edgewidth{edge}(diam1,diam1;l1,l2...)
% Here diam1 = measured diameter,
% l1 : contour length of the pt where the diameter was measured. 
% Manual entry doesn't store l1, l2 etc; stores default as 0
% Manual entry would store NT.edgewidth{edge}(diam1,diam2;0,0...) etc
% l1 and l2 are calculated from t_a * NT.edgelens{edge}
% where t_a is obtained from function distance2curve.m
%% plot img and network
fig = figure;
ax = axes(fig);


img = imread(imgfilename);
imshow(img)
hold on
options = struct();
options.edgeplotopt = {'LineWidth',4};
options.edgecolor = [0 1 0];
NT.plotNetwork(options);

set(gca,'YDir','reverse') %flip image to same orientation as img
hold on

for ec = 1:NT.nedge
    % show the edge widths that are already set on the network
    if ~(isempty(NT.edgewidth{ec}))
        for wc = 1:size(NT.edgewidth{ec},1)
            % intersection point
            pt = interp1(NT.cumedgelen{ec},NT.edgepath{ec},NT.edgewidth{ec}(wc,2));
            plot(pt(1),pt(2),'m.', 'MarkerSize',20)
        end
    end
end
hold off
%% control button

exitButton = uicontrol(gcf,'Style', 'pushbutton', ...
    'String', 'Exit drawline', ...
    'Position', [300 10 30 30], 'Callback',@setExitVal);

    function setExitVal(hObject,event)
        set(hObject, 'UserData', -1);
        hObject.UserData
        guidata(hObject);
    end
%%

answer = questdlg('Preferred way to enter diameters?', ...
    'Input diameters','Manual number entry','Measure Tool',...
    'Cancel', 'Cancel');
% Handle response
switch answer
    case 'Manual number entry'
        % manual entry subroutine
        measureTool = 0;
    case 'Measure Tool'
        % Measure using drawlines()
        measureTool = 1;        
    case 'Cancel'
        disp('Width entry subroutine ended')
        return
end
%%
if (measureTool)
    %%
    % Measure using drawlines()
    uiwait(msgbox('Draw lines across branch to measure diameter. Press any key to draw new line.', ...
        '','modal'));
    uiwait(msgbox('Click exit button or Esc to stop line draw', ...
        '','modal'));
    
    % Get all width inputs
    n = 1;
    % 27 : Esc
   
    exitVal = 1; 
    exitButton.UserData = 1;
    while (exitVal > 0) 
        n
        guidata(exitButton);
        exitVal = get(exitButton, 'UserData');
        if (exitVal == -1) 
            disp('Exit button terminated diameter entry') 
            break
        end
        
        disp('drawing line')
        lHandle(n) = drawline();

        
        k = waitforbuttonpress;
        if (k == 0) %if mouse key is detected, only exit loop when keyboard is pressed
            while (k == 0)
                %disp('clicked')
                guidata(exitButton);
                exitVal = get(exitButton,'UserData');
                if (exitVal == -1)
                    break
                else
                    k = waitforbuttonpress;
                end
                    
            end
            if (exitVal == -1)
                break
            end
        else
            value = double(get(gcf,'CurrentCharacter'));
        
            if (value == 27)
                disp('Esc terminated diameter entry')
                break
            end
        end
        n = n+1;
        %drawnow();
        pause(0.5);
        
        guidata(exitButton);
        exitVal = get(exitButton,'UserData');
        

    end
    %uiwait(msgbox('Hit enter when done adjusting', ...
    %    '','modal'));
    input('Hit enter when done');
    % remove deleted handles
    saveLines = lHandle(isvalid(lHandle));

    %nVals = zeros(NT.nedge,1);
    for n = 1:length(saveLines)
        % Assign to those edges whose paths cross lines
        % check which point on edgepath is it closest to
        % shortlist edgepaths first from nodepos
        coord = (saveLines(n).Position(1,:) + saveLines(n).Position(2,:))/2;
        [~, node] = min(sum((NT.nodepos - coord).^2,2));
        %node
        % check which edges connect to this node
        posEdges = 0;
        posEdges = [find(NT.edgenodes(:,2) == node); find(NT.edgenodes(:,1) == node)];
        checkpaths = NT.edgepath(posEdges);
        allPathPts = vertcat(checkpaths{:});
        pathLens = cell2mat(cellfun(@length,checkpaths,'uni',false));
        % get index where nearest path pt is located
        [~, pathind] = min(sum((allPathPts - coord).^2,2));
        edge = posEdges(find(cumsum(pathLens) >= pathind,1));
        %nVals(edge) = nVals(edge)+1;
        % get distance data
        edgewidth = sqrt(sum((saveLines(n).Position(2,:) - saveLines(n).Position(1,:)).^2));
        %NT.edgewidth{edge,nVals(edge)} = edgewidth;
        
        % interpolate coord along contour
        [~,~,t_a] = distance2curve(NT.edgepath{edge},coord);
        
        if (isempty(NT.edgewidth{edge}))
            NT.edgewidth{edge}(1,1) = edgewidth;
            % store contour length
            NT.edgewidth{edge}(1,2) = t_a * NT.edgelens(edge);
        else
            NT.edgewidth{edge}(end+1,1) = edgewidth;
            NT.edgewidth{edge}(end,2) = t_a * NT.edgelens(edge);
        end
        
        
        % save edgewidth
        saveedgewidth = NT.edgewidth;
        %%
    end   
    
    
else
    uiwait(msgbox('Click on edge to add two diameter values. Red = current active edge, Blue = edge already set manually', ...
        '','modal'));
    uiwait(msgbox('Click exit button or press Esc to leave program', ...
        '','modal'));
    exitVal = 1;
    button = 1;
    exitButton.UserData = 1;
    while (exitVal > 0)   % read ginputs until a mouse right-button occurs       
        [x,y,button] = ginput(1);
        exitVal = get(exitButton, 'UserData');
        if ((button == 27) | (exitVal < 0))
            disp('Exiting manual diameter entry program')
            uiwait(msgbox('Exiting manual diameter entry program', ...
                '','modal'));
            break
        end
        % check which point on edgepath is it closest to
        % shortlist edgepaths first from nodepos
        coord = [x,y];
        [~, node] = min(sum((NT.nodepos - coord).^2,2));
        %node
        % check which edges connect to this node
        posEdges = 0;
        posEdges = [find(NT.edgenodes(:,2) == node); find(NT.edgenodes(:,1) == node)];
        checkpaths = NT.edgepath(posEdges);
        allPathPts = vertcat(checkpaths{:});
        pathLens = cell2mat(cellfun(@length,checkpaths,'uni',false));
        % get index where nearest path pt is located
        [~, pathind] = min(sum((allPathPts - coord).^2,2));
        edge = posEdges(find(cumsum(pathLens) >= pathind,1));
        
        % highlight which edge is detected ultimately
        sn = plot(NT.edgepath{edge}(:,1),NT.edgepath{edge}(:,2),'color','r','Linewidth',2);
        prompt = {'Diameter1 after branching from mother'; 'Diameter2 before branching to daughter branches'};
        % diameter 1 : diameter after branching from mother,
        % diameter 2: diameter before branching to daughter branches
        dlgtitle = 'Enter diameter for red edge';
        dims = [1 40];
        if (isempty(NT.edgewidth{edge}))
            definput = {'',''};
        elseif (length(NT.edgewidth{edge}(:,1)) == 1)
            definput = {num2str(NT.edgewidth{edge}(1,1)),''};
        else
            definput = {num2str(NT.edgewidth{edge}(1,1)),num2str(NT.edgewidth{edge}(2,1))};
        end
        
        
        userInp = inputdlg(prompt,dlgtitle,dims,definput);
        edgewidth = cellfun(@str2num,userInp);
        delete(sn)
        plot(NT.edgepath{edge}(:,1),NT.edgepath{edge}(:,2),'color','b','Linewidth',2)
        % reset pre-existing edgeval
        NT.edgewidth{edge}(1,1) = edgewidth(1);
        NT.edgewidth{edge}(2,1) = edgewidth(2);
        saveedgewidth = NT.edgewidth;
        
        
        
        
    end
    saveedgewidth = NT.edgewidth;
end

end


