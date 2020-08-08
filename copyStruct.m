function opt2 = copyStruct(opt1,opt2,varargin)
% copy fields from structure opt1 to structure opt2
% if addnew=1, then add additional fields into opt2 from opt1
% exclude is a cell list of fields to exclude


addnew = 0;
exclude={};

for vc = 1:2:length(varargin)
    switch (varargin{vc})
        case('addnew')
            addnew = varargin{vc+1};
        case('exclude')
            exclude = varargin{vc+1};
    end
end

inputopt = fieldnames(opt1);
for c = 1:length(inputopt)
    s = inputopt(c); s=s{1};
    if ismember(s,exclude)
        continue
    end
    if (addnew || isfield(opt2,s))
        opt2.(s) = opt1.(s);
    end
end


end