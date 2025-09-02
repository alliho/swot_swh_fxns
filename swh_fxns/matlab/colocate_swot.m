function [obssw] = colocate_swot(SWOT, t, lon, lat, varargin)


%% set default variables [for co-location]
interpopt = 0; ri = 3;
Nt = length(t);
lon0 = nanmean(lon); lat0 = nanmean(lat);


%% parse variable inputs
p = inputParser; p.KeepUnmatched = true;
addParameter(p, 'ri', ri);
addParameter(p, 'interp', interpopt);
addParameter(p, 'prcnan', 0.5);


parse(p, varargin{:});
ri = p.Results.ri;
interpopt = p.Results.interp;
prcnan = p.Results.prcnan;

%% define fields - only ones that already exist
fldnms = fieldnames(SWOT);
dims = cellfun(@(fld) size(SWOT(1).(fld)) , fldnms, 'Un', 0);
ii = find(~cellfun(@(dim) sum(dim<=2), dims));
fldnms = fldnms(ii)';

fldnms = [{'time'}, fldnms];
fldnmsinit = fldnms;
fldnmsinit{2,1} = {};
fldnmsinit(2,:) = {NaN(size(SWOT))};

%% pull swot
clear obssw
obssw = struct(fldnmsinit{:});
for si=1:length(SWOT)
    swot = SWOT(si); if isempty(swot.latitude) | allnan(swot.time); continue; end
    X = swot.longitude; Y = swot.latitude; Nx = size(X,1); Ny = size(X,2);
    if Nt>1
        [dt, ti] = min(abs(t - nanmean(swot.time))); dt_o = mode(diff(t)); 
        if dt > dt_o/2; continue; end %skip if no data within a half hour of SWOT passover
        lon0 = lon(ti); lat0 = lat(ti);
    end

    ii = find(isin(X(1,:)-360, lon0 + [-1 1].*5) & isin(Y(1,:), lat0 + [-1 1].*5)); if isempty(ii); continue; end
    dist = distance(lat0, lon0 , Y(:,ii),X(:,ii) - 360, referenceEllipsoid('wgs84'))./1000;
    
    minValue = min(min(dist)); % just picking closest pixel to coordinates

    if minValue > 10; continue; end % if no close pixel then skip 


    [xi, yi] = find(dist==minValue); yi = ii(yi); 
    di = max([1 floor(ri/2)]); % radius to compute interpolation around 
    xx = xi-di:xi+di; yy = yi-di:yi+di; xx = xx(xx>0 & xx<=Nx); yy = yy(yy>0 & yy<=Ny);
 
    for fi=1:length(fldnms)
        fldnm = fldnms{fi};
        if strcmp(fldnm, 'time')
            obssw.(fldnm)(si) = swot.(fldnm)(yi); % ASSIGN
            continue;
        end

        var = swot.(fldnm);
        if contains(fldnm, 'ssha') & contains(fldnm, 'proc')
            var = swot.(fldnm) + swot.height_cor_xover;
        end


        if interpopt  
            try
                Xs = X(xx,yy) - 360; Ys = Y(xx,yy); Zs = var(xx,yy); 
            catch
                error('didntwork')
            end
            x0 = lon0; y0 = lat0;
            nn = find(~isnan(Zs(:))); 
            prc_nan = 1 - length(nn)./length(Zs(:));
            % obssw.([fldnm '_interp_prc_nan'])(si) = prc_nan; 
            if prc_nan>prcnan | length(nn)<2
                z0 = NaN;
            else
                F = scatteredInterpolant(Xs(nn),Ys(nn),Zs(nn));
                z0 = F(x0,y0);
            end
            obssw.(fldnm)(si) = z0; % ASSIGN        
        else
            obssw.(fldnm)(si) = var(xi,yi); % ASSIGN
        end

    end
end
   

    
end