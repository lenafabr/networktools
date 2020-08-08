function [nodepos,edgenodes,edgevals,nodelabels] = loadnetworkstruct(fname,options)
% load network structure from file fname

% load several structures at once?
opt.loadseveral=0;
% if loading several structures, use this word to note start of new
% structure
opt.structword = 'structure';

if (exist('options','var'))
    opt = copyStruct(options,opt);
end

%%
fid = fopen(fname);

tline = fgetl(fid);
clear nodepos edgenodes
struct = 0;

while (ischar(tline))
    
    checkstruct = strfind(tline,opt.structword);
    if (~isempty(checkstruct)) % new structure detected
        if (struct>0) % save already-found structure
            nodeposall{struct} = nodepos;
            edgenodesall{struct} = edgenodes;
        end
        % move on to next structure
        struct = struct + 1;
        if (mod(struct,10)==0); disp(struct); end
        if (~opt.loadseveral & struct > 1)
            break
        end         
    end
    
    % ignore comments and blank lines
    if (length(tline)==0 || tline(1)=='#'); tline = fgetl(fid); continue; end
    
    %disp(tline);
    
    % split line by whitespace
    words= strsplit(tline);
    
    if(strcmpi(words(1),'NODE'))
        % read in node position       
        %dim = nnz(cellfun(@(i) length(i), words(2:end))>0)-1;  
        %% matlab-compatible exponential mode
        for wc = 2:length(words)
            words{wc} = strrep(words{wc},'D','E');
        end
        nums= cellfun(@(i) str2double(i), words(2:end));
        nums = nums(~isnan(nums));
        dim = length(nums)-1;
        nodepos(nums(1),:) = nums(2:end);      
        if (length(words)>=dim+3)
            nodelabels{nums(1)} = words{dim+3};
        else
            nodelabels{nums(1)} = '';
        end
    elseif(strcmpi(words(1),'EDGE'))
        nums= cellfun(@(i) str2num(i), words(2:end));
        edgenodes(nums(1),:) = nums(2:3);
        edgevals(nums(1),:) = nums(4:end);
    end
    tline = fgetl(fid);
    
end

fclose(fid);

if (opt.loadseveral)
    nodepos = nodeposall;
    edgenodes = edgenodesall;
else
end

end