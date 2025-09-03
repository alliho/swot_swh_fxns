function [] = format_mapticks(ha, varargin)

if ~isempty(varargin)
    axstr = varargin{1}; 
else
    axstr = 'xy';
end

if length(axstr)==1
    for hi=1:length(ha)
        lbs = get(ha(hi), [upper(axstr) 'TickLabels']);
        lbs = cellfun(@(x) [x '^o'], lbs, 'Un',0);
        set(ha(hi), [upper(axstr) 'TickLabels'],lbs);
    end
elseif strcmp(axstr, 'xy')
    axstr = 'x';
    for hi=1:length(ha)
        lbs = get(ha(hi), [upper(axstr) 'TickLabels']);
        lbs = cellfun(@(x) [x '^o'], lbs, 'Un',0);
        set(ha(hi), [upper(axstr) 'TickLabels'],lbs);
    end
    axstr = 'y';
    for hi=1:length(ha)
        lbs = get(ha(hi), [upper(axstr) 'TickLabels']);
        lbs = cellfun(@(x) [x '^o'], lbs, 'Un',0);
        set(ha(hi), [upper(axstr) 'TickLabels'],lbs);
    end
end

end