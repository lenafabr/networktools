function setPlotDataTips(h,tipval)
% for plotted curves, set up the datatip info to show a particular value
% h = handle to the plotted line
% tipval = what value to show along this curve

dttemplate = h.DataTipTemplate;
dttemplate.FontSize=6;
dttemplate.DataTipRows(1).Value = tipval*ones(size(h.XData));
dttemplate.DataTipRows(1).Label = '';
dttemplate.DataTipRows(2:end) = [];
end