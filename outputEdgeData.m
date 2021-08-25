function outputEdgeData(NT,outputfile)
% this function writes out the edge lengths and widths to a file
% widths (along a *single* edge) are ordered from upstream to downstream measurement 

fid = fopen(outputfile,'w');

writeline ='Edge ID, length, width (ordered from upstream to downstream)';
disp(writeline);
fprintf(fid,[writeline '\n']);

for ec = 1:NT.nedge
    if (size(NT.edgewidth{ec},2)>1)
        [~,ind] = sort(NT.edgewidth{ec}(:,2));        
        widths = NT.edgewidth{ec}(ind,1);
    elseif(~isempty(NT.edgewidth{ec}))
        widths = NT.edgewidth{ec}(:,1);
    else
        widths=[];
    end
    widthstring = sprintf('%f, ',widths');
    writeline = sprintf('%d %f %s',ec, NT.edgelens(ec), widthstring(1:end-2));
    disp(writeline);
    fprintf(fid,[writeline '\n']);
end
fclose(fid);

end