function [wb] = load_any_nc(fname, varargin)
% -------------------------------------------------------------------------
% LOAD_ANY_NC  Loads all variables of netcdf into struct 
% Alli Ho 2019
% -------------------------------------------------------------------------
%% GET INFO
info = ncinfo(fname);
if isempty(info.Variables) & ~isempty(info.Groups); 
    groups = 1; 
else
    groups = 0; 
end
%% OPTIONS 
limsargin = {};
dims = {};

p = inputParser; 
addParameter(p, 'groups', groups);
addParameter(p, 'group', []);
addParameter(p, 'Dimensions', dims);
addParameter(p, 'Range', limsargin);
parse(p, varargin{:});
groups = p.Results.groups;
grpnm = p.Results.group;
limsargin = p.Results.Range;
dims = p.Results.Dimensions;

if ~(~isempty(dims) & ~isempty(limsargin))
    limsargin = {};
    dims = {};
end



%% LOAD

if ~groups
    % variable list

    paramlist = strings(1,length(info.Variables)); 
    for k=1:length(info.Variables)
        paramlist(k) = info.Variables(k).Name;
    end
elseif groups
    grps = {info.Groups.Name};
    if ~isempty(grpnm)
        gg = find(cellfun(@(x) strcmp(grpnm, x), grps));
        grps = grps(gg); 
    end
    Ng = length(grps);

    paramlists = [];
    for i=1:Ng
        info = ncinfo(fname, grps{i});
        paramlist = strings(1,length(info(1).Variables)); 
        for k=1:length(info.Variables)
            paramlist(k) = ['/' grps{i} '/' info.Variables(k).Name];
        end
        paramlists = [paramlists paramlist];
    end
    paramlist = paramlists;
end




% load variables into structure
wb = load_nc_from_paramlist([fname], paramlist, limsargin{:}, 'Dimensions', dims);  % function created to load netCDF into structs


% add time unit if variable
if sum(contains(paramlist, 'time')) & ~groups
    pi = find(strcmp(paramlist, 'time')); 
    if isempty(pi); pi = find(contains(paramlist, 'time')); end
    if length(pi)>1; pi = pi(1); end
    ui = find(contains({info.Variables(pi).Attributes(:).Name}, 'units'));
    if isempty(ui); return; end
    unitstr = info.Variables(pi).Attributes(ui).Value;
    wb.timeinfo.longname = unitstr; 
    if contains(unitstr, 'day'); wb.timeinfo.unit = 'day'; 
    elseif contains(unitstr, 'second'); wb.timeinfo.unit = 'second'; 
    elseif contains(unitstr, 'hour'); wb.timeinfo.unit = 'hour'; end
    
    unitstr = strsplit(unitstr, ' since '); 
    if contains(unitstr{2},'UTC'); unitstr{2} = strrep(unitstr{2}, 'UTC', ''); end
    if contains(unitstr{2},'GMT'); unitstr{2} = strrep(unitstr{2}, 'GMT', ''); end
    if contains(unitstr{2},'T'); t0 = datenum(unitstr{2}, 'yyyy-mm-ddTHH:MM:SS'); 
    else; t0 = datenum(unitstr{2}); end
    
    wb.timeinfo.t0 = t0;
end


end

function [ outlist ] = load_nc_from_paramlist( filename, paramlist, varargin)
    % -------------------------------------------------------------------------
    % LOAD_NC_FROM_PARAMLIST  Loads netcdf into struct given parameter list
    % Alli Ho 2023 based on LOADNCDF code from 2019
    % varargin has to be inputs to ncread (start, count, stride)
    % -------------------------------------------------------------------------
    subsetopt = 0;
    if ~isempty(varargin)
        if length(varargin) > 2
            subsetopt = 1;
        end
    end
    
    if subsetopt
        if sum(strcmp(varargin, 'Dimensions'))
            di = find(strcmp(varargin, 'Dimensions')); 
            dims = varargin{di+1};

            % get dimension names
            dinfo = cellfun(@(param) ncinfo(filename, param).Dimensions, paramlist, 'Un', 0);
            ii = find(~cellfun(@isempty, dinfo));
            dinfo = cellfun(@(param) [ncinfo(filename, param).Dimensions.Length], paramlist(ii), 'Un', 0);
            dd = find(cellfun(@length, dinfo)==length(dims));
            dd = dd(find(cellfun(@(x) sum(x==dims), dinfo(dd))));
            pi = ii(dd); pi = pi(1);
            dimnames = {ncinfo(filename, paramlist{pi}).Dimensions.Name};
        end
        [start, count, stride] = varargin{1:3};
    end



    %%% check if only one paramlist group
    paramheader = cellfun(@(x) strsplit(x, '/'), cellstr(paramlist), 'Un', 0);
    if ~all(cellfun(@length, paramheader)>1); paramheader = NaN;
    else
        paramheader = cellfun(@(x) x{2}, paramheader, 'Un', 0);
        if all(cellfun(@(x) strcmp(x, paramheader{1}), paramheader));
            paramheader = paramheader{1}; 
        else
            paramheader = NaN;
        end
    end

    for n=1:length(paramlist)
        thisparam = char(paramlist(n));
        
        fieldname = thisparam; 
        fieldname = strrep(fieldname, '/', '_');
        fieldname = strrep(fieldname, '.', '_');
        if ~isnan(paramheader); 
            fieldname = strrep(fieldname, paramheader, '');
        end

        while ~isletter(fieldname(1))
            fieldname = fieldname(2:end);
        end
        
        info = ncinfo(filename, thisparam); 
        if subsetopt & ~isempty(info.Dimensions)
            thisdims = [info.Dimensions.Length];
            % thisdims = [info.Size];
            thisdimnames = {info.Dimensions.Name};

            % [~, di_d] = intersect(thisdims, dims);
            di_d = find(sum(thisdims'==dims,1));
            di_n = find(ismember(dimnames, thisdimnames));
            di = intersect(di_d, di_n);
            % di = find(thisdims==dims & ismember(dimnames, thisdimnames)); %di = 1:di;
            limsargin = {start(di), count(di),  stride(di)};
            try
                S1.(fieldname) = ncread(filename, thisparam, limsargin{:});
            catch
                % disp(['ERROR : range limit not applicable to ' thisparam]); 
                % disp(['        dimensions [' num2str(thisdims) ']']); 
            end
        else
            S1.(fieldname) = ncread(filename, thisparam);
        end
        
    end
    outlist = S1;
end

