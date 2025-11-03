function [SWOT] = load_swot(fpath, varargin)



%% SET PARAMETERS
%%% set file heirarchies
fproc_heirarchy = {'PIC2', 'PGC0', 'PIC0', 'PID0'};
add_expert_fields = 0; 
fstr = '*SWOT*.nc';
verbose = 1;
%% variable inputs


if ~isempty(varargin)
    if length(varargin)==2 & sum(cellfun(@isnumeric,varargin(1:2)))==2
        [lonlims, latlims] = varargin{1:2};
    elseif length(varargin)>=3 
        if sum(cellfun(@isnumeric,varargin(1:3)))==3
            [lonlims, latlims, daterange] = varargin{1:3};
            varargin = varargin(4:end);
        elseif sum(cellfun(@isnumeric,varargin(1:3)))==2
            [lonlims, latlims] = varargin{1:2};
            varargin = varargin(3:end);
        end
    end
end

p = inputParser; p.KeepUnmatched = true;
addParameter(p, 'search', fstr);
addParameter(p, 'pass', NaN);
addParameter(p, 'cycle', NaN);
addParameter(p, 'verbose', 1);
addParameter(p, 'omit', {});
addParameter(p, 'include', {});
addParameter(p, 'add', '');
addParameter(p, 'combine', 0);
parse(p, varargin{:});
fstr = p.Results.search;
choosepass = p.Results.pass;
choosecycle = p.Results.cycle;
verbose = p.Results.verbose;
omitflds = p.Results.omit;
incflds = p.Results.include;
addopt = p.Results.add;
combineswaths = p.Results.combine;

if ~isempty(addopt)
    if contains(addopt, 'expert')
        add_expert_fields = 1;
    end
end


if isnumeric(choosecycle)
    if ~isnan(choosecycle); 
    choosecycle = ['000' num2str(choosecycle)]; choosecycle = choosecycle(end-2:end);
    end
end
if ~isnumeric(verbose); if strcmp(verbose, 'off'); verbose = 0; elseif strcmp(verbose, 'on'); verbose = 1; end; end

if ~exist('latlims')
    latlims = [-80 80]; lonlims = [-360 360];
end

swottype = 'Expert';
if contains(fstr, 'Unsmoothed'); swottype = 'Unsmoothed'; end

% if ~isempty(incflds)
% 
% end


if ~strcmp(swottype, 'Unsmoothed'); combineswaths = 0; end

% if ~isempty(varargin)
%     if length(varargin)==2 & sum(cellfun(@isnumeric,varargin(1:2)))==2
%         [lonlims, latlims] = varargin{1:2};
%     elseif length(varargin)>=3 & sum(cellfun(@isnumeric,varargin(1:3)))==3
%         [lonlims, latlims, daterange] = varargin{1:3};
%     end
% 
% 
% 
%     vi = find(strcmp(varargin, 'search'));
%     if ~isempty(vi)
%         vi = vi(end);
%         fstr = varargin{vi+1};
%     end
%     vi = find(strcmp(varargin, 'pass'));
%     if ~isempty(vi)
%         choosepass = varargin{vi+1};
%     end
%     vi = find(strcmp(varargin, 'cycle'));
%     if ~isempty(vi)
%         choosecycle = varargin{vi+1};
%         if isnumeric(choosecycle);
%             choosecycle = ['000' num2str(choosecycle)]; choosecycle = choosecycle(end-2:end);
%         end
%     end
%     vi = find(strcmp(varargin, 'verbose'));
%     if ~isempty(vi)
%         verbose = varargin{vi+1};
%         if ~isnumeric(verbose); if strcmp(verbose, 'off'); verbose = 0; elseif strcmp(verbose, 'on'); verbose = 1; end; end
%     end
% end
% if ~exist('latlims')
%     latlims = [-80 80]; lonlims = [-360 360];
% end

%% get file names
fnames = dir([fpath fstr]);
ftimes = cellfun(@(x) x{8}, cellfun(@(x) strsplit(x, '_'), {fnames.name}, 'Un', 0), 'Un', 0);
ftimes = cellfun(@(x) datenum(x, 'yyyymmddTHHMMSS'), ftimes);
fcycle = cellfun(@(x) str2num(x{6}), cellfun(@(x) strsplit(x, '_'), {fnames.name}, 'Un', 0));
fpass = cellfun(@(x) str2num(x{7}), cellfun(@(x) strsplit(x, '_'), {fnames.name}, 'Un', 0));
fproc = cellfun(@(x) x{10}, cellfun(@(x) strsplit(x, '_'), {fnames.name}, 'Un', 0), 'Un', 0);
fext = cellfun(@(x) x{1}(end-6:end), cellfun(@(x) strsplit(x, '.'), {fnames.name}, 'Un', 0), 'Un', 0);

%% SUBET by time
if ~exist('daterange') | ~isnumeric(daterange)
    daterange = minmax(ftimes) + [-1 1]; 
    if allnan(daterange)
        daterange = nanmean(ftimes) + [-1 1]; 
    end
end

tt = find(isin(ftimes, daterange));
if ~isnan(choosepass)
    tt = intersect(tt, find(fpass==str2num(choosepass)));
end
if ~isnan(choosecycle)
    tt = intersect(tt, find(fcycle==str2num(choosecycle)));
end
if isempty(tt)
    error('no files available')
end
A = {fnames, ftimes, fcycle, fpass, fproc, fext};
A = cellfun(@(x) x(tt), A, 'Un', 0); 
[fnames, ftimes, fcycle, fpass, fproc, fext] = A{:};

%%
if sum(cellfun(@(x) contains(x, 'Expert'), {fnames.name}))==length({fnames.name})
    swottype = 'Expert';
elseif sum(cellfun(@(x) contains(x, 'Unsmoothed'), {fnames.name}))==length({fnames.name})
    swottype = 'Unsmoothed'; 
else
    error('no or mixed swottype');
end

%% load SWOT
utimes = unique(ftimes); 
N = length(utimes); 
si = 1; 
ti = 1; 
clear SWOT
while ti <= N
    t0 = utimes(ti); 

    dispv(1, ['time: ' datestr(t0)])
    dispv(verbose, ['      swot index: ' num2str(si)])
    dispv(1, ['      time index: ' num2str(ti) '/' num2str(N)])

    ff = find(ftimes==t0); FF = ff; 
    Fexts = cellfun(@(x) str2num(x(end-4:end-3)), {fnames(FF).name}); 
    Fprocs = fproc(FF);
    [~, ~, miproc] = intersect(Fprocs, fproc_heirarchy);
    mpi = find(strcmp(Fprocs, fproc_heirarchy(min(miproc))));
    [~, mai] = max(Fexts(mpi)); mi = mpi(mai); 
    fi = FF(mi);

    if length(FF)>1 % say what happened if you had to choose a file
        dispv(verbose, ['      extensions:']); cellfun(@(x) dispv(verbose, ['      - ' x]), fext(FF))
        dispv(verbose, ['      ---> chose file with tail: ' fnames(fi).name(end-9:end-3)])
    end
    dispv(verbose, ['      file index: ' num2str(fi) '/' num2str(length(fnames))])
    tic
    


    %%% TEST TO SEE IF WITHIN LIMITS and CHOOSE RANGE TO LOAD
    limsargin = {};
    try
        if strcmp(swottype, 'Expert')
            lon = ncread([fpath fnames(fi).name], 'longitude')-360;
            lat = ncread([fpath fnames(fi).name], 'latitude');
            info = ncinfo([fpath fnames(fi).name], 'latitude'); dims = [info.Dimensions.Length]; 
        elseif strcmp(swottype, 'Unsmoothed')
            lon = ncread([fpath fnames(fi).name], '/right/longitude')-360;
            lat = ncread([fpath fnames(fi).name], '/right/latitude');
            info = ncinfo([fpath fnames(fi).name], '/right/latitude'); dims = [info.Dimensions.Length]; 
        end
        inbounds = (isin(lon, lonlims + [-1 1].*0.5)  & isin(lat, latlims + [-1 1].*0.5));
        yy = find(sum(inbounds,1)); xx = find(sum(inbounds,2));
        start = [min(xx) min(yy)]; count = [range(xx) range(yy)]+1; stride = [1 1];

        % if isempty(xx) | isempty(yy)
        % 
        % 
        % end

        % don't limit xx (cross track range)
        start(1) = 1; count(1) = dims(1);

        strcntstrd = {start; count; stride};
        
        if sum(count<=1)
             ti = ti + 1; dispv(verbose, ['      ERROR: outside buoy limits *** *** ***']); 
             continue;
        end
        limsargin = {'Dimensions', dims, 'Range', strcntstrd};
        dispv(1, ['      --> extracting: [' num2str(count(1)) ',' num2str(count(2))  '] out of [' num2str(dims(1)) ',' num2str(dims(2)) ']'])
        
    catch
        % ti = ti + 1; 
        dispv(verbose, ['      ISSUE:  assigning limits *** *** ***']);
    end


    %%% NOW FINALLY LOADING SWOT! 
    try
        swot = load_any_nc([fpath fnames(fi).name], limsargin{:}, 'omit', omitflds); 
    catch
        ti = ti + 1; 
        dispv(verbose, ['      ERROR: data read error *** *** ***']);
        continue; 
    end

    if strcmp(swottype, 'Expert')
        % seconds since 2000-01-01 00:00:00.0
        swot.time = swot.time./(3600*24) + datenum(2000,1,1);
        % filter
        Ny = size(swot.latitude,2);
        yy = find( isin(swot.latitude(1,:) ,latlims + [-1 1].*0.5) & isin(swot.longitude(1,:), (lonlims  + [-1 1].*0.5) + 360) | isin(swot.latitude(end,:) ,latlims + [-1 1].*0.5) & isin(swot.longitude(end,:), (lonlims  + [-1 1].*0.5) + 360) ); 
        
        if length(yy)*2 < 60
            ti = ti + 1; 
            dispv(verbose, ['      ERROR: not enough observations *** *** ***']);
            continue
        end
    elseif strcmp(swottype, 'Unsmoothed')
        fldnms = fieldnames(swot);
        ii = find(cellfun(@(x) strcmp(x(end-3:end), 'time'), fldnms)); 
        tfldnms = fldnms(ii); 
        for i=1:length(tfldnms)
            swot.(tfldnms{i}) = swot.(tfldnms{i})./(3600*24) + datenum(2000,1,1);
        end
        ii = find(cellfun(@(x) strcmp(x(end-7:end), 'latitude'), fldnms)); fldlat = (fldnms{ii(1)}); 
        ii = find(cellfun(@(x) strcmp(x(end-8:end), 'longitude'), fldnms)); fldlon = (fldnms{ii(1)}); 
        
        Ny = size(swot.(fldlat),2); Nx = size(swot.(fldlat),1);
        xi = floor(Nx/2);
        yy = find(isin(swot.(fldlat)(xi,:) ,latlims + [-1 1].*1) & isin(swot.(fldlon)(xi,:), (lonlims  + [-1 1].*1) + 360)); 
        swot.time = nanmean(cell2mat(cellfun(@(fld) swot.(fld)', tfldnms, 'Un', 0)));
    end
    
    fldnms = fieldnames(swot);
    for ei = 1:length(fldnms)
        fldnm = fldnms{ei};
        tmp = swot.(fldnm);

        if size(tmp,1) == Ny
            swot.(fldnm) = tmp(yy,:);
        elseif size(tmp,2) == Ny
            swot.(fldnm) = tmp(:,yy); 
        end

    end


   




    if isempty(swot.time); ti = ti + 1; dispv(verbose, ['      ERROR: no data *** *** ***']); continue; end
    swot.fname = fnames(fi).name;
    swot.fpath = fpath;
    try
        [swot.t0, swot.cycle, swot.pass, swot.processing, swot.orientation, swot.angle] = get_swot_info(swot, swot.fname);
    end

    %%% combining left and right swaths
    if strcmp(swottype, 'Unsmoothed') & combineswaths
        swot = combine_leftright(swot);
        swot = rmfield(swot, {'left_polarization_karin','left_time_tai', 'right_polarization_karin', 'right_time_tai'});
    end

    SWOT(si) = swot; 
    dispv(verbose, ['      ASSIGNED to SWOT(' num2str(si) ')']);
    si = si + 1; ti = ti + 1; % saved
    toc
    
end


if ~exist('SWOT'); error('NO DATA AVAILABLE'); end

ii = find(~cellfun(@isempty, {SWOT(:).time}));
SWOT = SWOT(ii);


%% adding expert fields to unsmoothed


if strcmp(swottype, 'Unsmoothed') & add_expert_fields
    disp('----- adding Expert fields to Unsmoothed ')
    %%% only field names:
    fldnms = {'height_cor_xover','mean_sea_surface_cnescls', 'solid_earth_tide', 'ocean_tide_fes', 'internal_tide_hret','pole_tide', 'dac', 'cross_track_distance'};

    fstr = strrep(fstr, 'Unsmoothed', 'Expert');
    fpath = strrep(fpath, 'unsmoothed/', '');
    [SWOT_ex] = load_swot(fpath, lonlims, latlims, daterange, varargin{:}, 'search', fstr, 'verbose',0, 'include', [{'latitude', 'longitude', 'time'} fldnms]);
    t0_ex = [SWOT_ex.t0];
    t0 = [SWOT.t0];
    [~, ii] = intersect(round(t0.*24)./24, round(t0_ex.*24)./24);
    disp('----- interpolating')
    
    if ~isempty(ii) 
        for si=ii'
            swot_un = SWOT(si);
            six = find(round(t0_ex.*24)./24==round(swot_un.t0.*24)./24);
            swot = SWOT_ex(six);
    
            %%% interpolate expert to unsmoothed
            if sum(strcmp(fieldnames(swot), 'right_') | strcmp(fieldnames(swot), 'left_')) > 0.5.*length(fieldnames(swot))
                grps = {'left_', 'right_'};
            else
                grps = {''};
            end

            X2 = swot.longitude;
            Y2 = swot.latitude;
            tic
            for gi=1:length(grps)
                grp = grps{gi};
                X1 = swot_un.([grp 'longitude']);
                Y1 = swot_un.([grp 'latitude']);
                tmp = ones(size(Y2));
                nn = find(isin(Y2(:), minmax(Y1) + [-1 1].*0.05));
                handshake = scatteredInterpolant(X2(nn), Y2(nn),tmp(nn));


                for fli=1:length(fldnms)
                    fldnm = fldnms{fli};
                    Z2 = swot.(fldnm);
                    
                    
                    % nn = find(~isnan(Z2(:)) & isin(Y2(:), minmax(Y1) + [-1 1].*0.05));
                    % nn = find(~isnan(Z2(:)));
                    % handshake.Values = 
                    % Z_F = handshak
                    Z_F = handshake;
                    Z_F.Values = Z2(nn);


                    % Z_F = scatteredInterpolant(X2(nn), Y2(nn),Z2(nn));
                    Z1 = Z_F(X1, Y1);


                    % [~, ti_mi] = min(abs(t(ti) - tbins));
                    % Z = swotcube(pai).(fldnm)(:,:,ti);
            
                    % [ matrix, nmatrix, xmid, ymid] = bin2D(xbins, ybins, X,Y,Z);
                    % Zb(1:end-1,1:end-1,ti_mi) = matrix;
            
                    % nn = find(~isnan(Z(:)) & isin(Z(:), prctile(Z(:), prclims)));
                    % [N,Xedges,Yedges,binX,binY] = histcounts2(X(nn),Y(nn),xbins,ybins);
                    % 
                    % 
                    % ii = find(binX~=0 & binY~=0); if isempty(ii); continue; end
                    % [nx ny] = size(N);
                    % if contains(fldnm, 'direction') | contains(fldnm, 'delta_currvswv')
                    %     Zg = accumarray([binX(ii),binY(ii)],Z(nn(ii)),[nx ny],@dirnanmean,NaN);
                    % else
                    %     Zg = accumarray([binX(ii),binY(ii)],Z(nn(ii)),[nx ny],@mean,NaN);
                    % end
                    % Zg(end+1,end+1) = NaN;
                    % Zg(end,:) = NaN;
                    % Zg(:,end) = NaN;
                    % n(1:size(N,1),1:size(N,2)) = N;
                    % 



                    SWOT(si).([grp fldnm]) = Z1;
                end
            end
            toc
    
        end
    end
    toc
end



end

function dispv(verbose, str)
    if verbose
        disp(str);
    end
end


