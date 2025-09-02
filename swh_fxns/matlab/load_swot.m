function [SWOT] = load_swot(fpath, varargin)



%% SET PARAMETERS
%%% set file heirarchies
fproc_heirarchy = {'PIC2', 'PGC0', 'PIC0', 'PID0'};
add_expert_fields = 0; 
fstr = '*SWOT*.nc';

%% variable inputs
if ~isempty(varargin)
    if length(varargin)==2 & sum(cellfun(@isnumeric,varargin))
        [lonlims, latlims] = varargin{1:2};
    elseif length(varargin)>=3 & sum(cellfun(@isnumeric,varargin))
        [lonlims, latlims, daterange] = varargin{1:3};
    end
    vi = find(strcmp(varargin, 'search'));
    if ~isempty(vi)
        vi = vi(end);
        fstr = varargin{vi+1};
    end
end
if ~exist('latlims')
    latlims = [-80 80]; lonlims = [-360 360];
end

swottype = 'Expert';
if contains(fstr, 'Unsmoothed'); swottype = 'Unsmoothed'; end

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

A = {fnames, ftimes, fcycle, fpass, fproc, fext};
tt = find(isin(ftimes, daterange));
A = cellfun(@(x) x(tt), A, 'Un', 0); 
[fnames, ftimes, fcycle, fpass, fproc, fext] = A{:};

%% load SWOT
utimes = unique(ftimes); 
N = length(utimes); 
si = 1; 
ti = 1; 
clear SWOT
while ti <= N
    t0 = utimes(ti); 

    disp(['time: ' datestr(t0)])
    disp(['      swot index: ' num2str(si)])
    disp(['      time index: ' num2str(ti) '/' num2str(N)])

    ff = find(ftimes==t0); FF = ff; 
    Fexts = cellfun(@(x) str2num(x(end-4:end-3)), {fnames(FF).name}); 
    Fprocs = fproc(FF);
    [~, ~, miproc] = intersect(Fprocs, fproc_heirarchy);
    mpi = find(strcmp(Fprocs, fproc_heirarchy(min(miproc))));
    [~, mai] = max(Fexts(mpi)); mi = mpi(mai); 
    fi = FF(mi);

    if length(FF)>1 % say what happened if you had to choose a file
        disp(['      extensions:']); cellfun(@(x) disp(['      - ' x]), fext(FF))
        disp(['      ---> chose file with tail: ' fnames(fi).name(end-9:end-3)])
    end
    disp(['      file index: ' num2str(fi) '/' num2str(length(fnames))])
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
        % don't limit xx (cross track range)
        start(1) = 1; count(1) = dims(1);

        strcntstrd = {start; count; stride};
        
        if sum(count<=1)
             ti = ti + 1; disp(['      ERROR: outside buoy limits *** *** ***']); 
             continue;
        end
        limsargin = {'Dimensions', dims, 'Range', strcntstrd};
        disp(['      --> extracting: [' num2str(count(1)) ',' num2str(count(2))  '] out of [' num2str(dims(1)) ',' num2str(dims(2)) ']'])
        
    catch
        % ti = ti + 1; 
        disp(['      ISSUE:  assigning limits *** *** ***']);
    end


    %%% NOW FINALLY LOADING SWOT! 
    try
        swot = load_any_nc([fpath fnames(fi).name], limsargin{:}); 
    catch
        ti = ti + 1; 
        disp(['      ERROR: data read error *** *** ***']);
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
            disp(['      ERROR: not enough observations *** *** ***']);
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





    if isempty(swot.time); ti = ti + 1; disp(['      ERROR: no data *** *** ***']); continue; end
    swot.fname = fnames(fi).name;
    try
        [swot.t0, swot.cycle, swot.pass, swot.processing, swot.orientation, swot.angle] = get_swot_info(swot, swot.fname);
    end

    SWOT(si) = swot; 
    disp(['      ASSIGNED to SWOT(' num2str(si) ')']);
    si = si + 1; ti = ti + 1; % saved
    toc
    
end



ii = find(~cellfun(@isempty, {SWOT(:).time}));
SWOT = SWOT(ii);



if strcmp(swottype, 'Unsmoothed') & add_expert_fields
    fstr = strrep(fstr, 'Unsmoothed', 'Expert');
    [SWOT_ex] = load_swot(fpath, varargin{:}, 'search', fstr);
    t0_ex = [SWOT_ex.t0];
    t0 = [SWOT.t0];
    [~, ii] = intersect(round(t0.*24)./24, round(t0_ex.*24)./24);
    if ~isempty(ii) 
        for si=ii'
            swot_un = SWOT(si);
            six = find(round(t0_ex.*24)./24==round(swot_un.t0.*24)./24);
            swot = SWOT_ex(six);
    
            %%% interpolate expert to unsmoothed
            grps = {'left', 'right'};
            fldnms = {'height_cor_xover','mean_sea_surface_cnescls', 'solid_earth_tide', 'ocean_tide_fes', 'internal_tide_hret','pole_tide', 'dac', 'cross_track_distance'};
            X2 = swot.longitude;
            Y2 = swot.latitude;
            for fli=1:length(fldnms)
                fldnm = fldnms{fli};
                Z2 = swot.(fldnm);
                for gi=1:length(grps)
                    grp = grps{gi};
                    X1 = swot_un.([grp '_longitude']);
                    Y1 = swot_un.([grp '_latitude']);
                    nn = find(~isnan(Z2(:)) & isin(Y2(:), minmax(Y1) + [-1 1].*0.2));
                    Z_F = scatteredInterpolant(X2(nn), Y2(nn),Z2(nn));
                    Z1 = Z_F(X1, Y1);
                    SWOT(si).([grp '_' fldnm]) = Z1;
                end
            end
    
        end
    end
end


end