base_path = '/Users/ajho/Documents/JPL/';
local_path = [base_path 'papers/swot_swh_calval/swot_swh_fxns/'];
addpath([local_path 'swh_fxns/matlab/'])
swot_fpath = [base_path '/data/SWOT/onedayrepeat/'];

%%% ADD DATA
load('/Users/ajho/Documents/myrepos/supportingdata/coastline_labeled.mat');
load('/Users/ajho/Documents/myrepos/supportingdata/nicejet.mat');

%%% OTHER DEPENDENCIES:
% m_map
% gitsio/code_universal/*/
    % isin.m
    % plot_mcoast.m
    % load_any_nc.m
    % setfigsize.m
    % drawbox.m
    % cleanLegend.m
    % quietaxes.m
    % textbypos.m
    % setaxes.m

% tight_subplot.m
% AddLetters2Plots.m


%% load

clear data
fpath = [local_path '/data/'];
tmp = load([fpath 'NDBC_x_SWOT_science.mat']); data(1) = tmp.data;
tmp = load([fpath 'GNSS_x_SWOT_calval.mat']); data(2) = tmp.data;
tmp = load([fpath 'NDBC_x_SWOT_calval.mat']); data(3) = tmp.data;
data(1).obstype = 'buoys';
data(2).obstype = 'moorings';
data(3).obstype = 'buoys';

%% define stats


getnn = @(x,y) find(~isnan(x) & ~isnan(y));
getnnn = @(x,y,z) find(~isnan(x) & ~isnan(y) & ~isnan(z));
biasfxn = @(x,y) nanmean((y-x));
rmsefxn = @(x,y) sqrt(nanmean((x-y).^2));
crmsefxn = @(x,y) sqrt(rmsefxn(x,y).^2 - biasfxn(x,y).^2);
normrmsefxn = @(x,y) sqrt(nanmean((x-y).^2))./nanmean(abs(x));
absrmsefxn = @(x,y) sqrt(nanmean((abs(x)-abs(y)).^2));
nfxn = @(x,y) sum(~isnan(x) & ~isnan(y));
ccfxn = @(x,y) mode(diag(flipud(corrcoef(x(getnn(x,y)),y(getnn(x,y))))));
sifxn = @(x,y) 1./nanmean(x(getnn(x,y))) .* sqrt( 1./length(x(getnn(x,y))) .* nansum( ( (y(getnn(x,y))-nanmean(y(getnn(x,y)))) - (x(getnn(x,y))-nanmean(x(getnn(x,y))))   ).^2 ) );;
linfxn = @(x,y) polyfit(x(getnn(x,y)), y(getnn(x,y)), 1);
getind = @(var,i) var(i); 
intcptfxn = @(x,y) getind(linfxn(x,y),2); 
slpfxn = @(x,y) getind(linfxn(x,y),1); 
slp0fxn = @(x,y) y(getnn(x,y))/x(getnn(x,y));

%% [save for paper] JUST SCATTER | buoys and moorings | with nadir? no model 
mdl = correct_swotswh('data', 'buoys');

dd = 1; thl = datenum('20240521T150311', 'yyyymmddTHHMMSS'); thl_lims = thl + [-1 1].*24/24; clims = [1.6 4.6];
thl = thl - 1;
dd = 2; thl = datenum(2023,4,19, 5, 34, 0); clims = [2.7 3.8]; thl_lims = thl + [-1 1].*1/24;


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
correctopt = 1; 
procname = 'P..0';
% procname = 'P...';
% procname = 'PIC2';
ogfldnms = {'swh_nadir_altimeter','swh_nadir_altimeter', 'swh_karin_uncorr_patch', 'swh_karin_uncorr_patch'};
% ogfldnms = {'swh_nadir_altimeter','swh_nadir_altimeter', 'swh_karin_uncorr_patch_unsmoothed', 'swh_karin_uncorr_patch_unsmoothed'};
hlopt = 1;
mdl = correct_swotswh();
mdl = correct_swotswh('data', 'moorings');
% mdl = correct_swotswh('data', 'buoys');


% -------------------------------------------------------------------------
dds = [dd dd];
if dd==2
    dlabs = {'GNSS buoy','GNSS buoy'};
else
    dlabs = {'coastal buoy','coastal buoy'};
end
fldnms = {'swh_karin_uncorr_patch', 'swh_model'};
% fldnms = {'swh_karin_uncorr_patch_unsmoothed', 'swh_model'};
% fldnms = {'swh_karin_proc_patch', 'swh_model'};
fldlbs = {'SWOT KaRIn', 'model'};

% -------------------------------------------------------------------------
figure(404); clf; 
setfigsize(gcf, [330   627])
% tight_subplot(Nh, Nw, [gap_h gap_w], [lower upper], [left right]) 
ha = tight_subplot(2, 1, [0.02 0.02], [0.08 0.015], [0.13 0.01]);

% -------------------------------------------------------------------------
setaxes(ha(1:end),'XLim', [0 7.5]);
setaxes(ha(1:end),'YLim', [0 7.5]);
% setaxes(ha(1:end),'XLim', [0 7.5]);
% setaxes(ha(1:end),'YLim', [0 7.5]);
setaxes(ha(:),'CLim', clims);
setaxes(ha(:),'Colormap', subsetcmap(nicejet, 20));
format_fig4print(ha)
setaxes(ha(:),'DataAspectRatio', [1 1 1]);
tx = [0:10];
setaxes(ha, 'XTick', tx)
setaxes(ha, 'YTick', tx)

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

% qcflds = {'swh_karin_proc'};



for hi=1:length(ha)
    axes(ha(hi)); hold on; 
    tx = [0:15];
    plot(tx,tx,'k--')
    % plot(tx,tx*1.05,'k:')
    % plot(tx,tx*0.95,'k:')
end



for hi=1:length(ha)
    disp(hi)
    dd = dds(hi);
    xlb = dlabs{hi};
    fldnm = fldnms{hi};
    ylb = fldlbs{hi};
    % -------------------------------------------------------------------------
    x = cell2mat(cellfun(@(x) x(:)', {data(dd).hs}, 'Un', 0));
    y = cell2mat(cellfun(@(x) x(:)', {data(dd).(fldnm)}, 'Un', 0)); y(y<0) = NaN;
    z = cell2mat(cellfun(@(x) x(:)', {data(dd).cross_track_distance}, 'Un', 0)); z = abs(z);

    t = cell2mat(cellfun(@(x) x(:)', {data(dd).t}, 'Un', 0));
    lat = cell2mat(cellfun(@(x) x(:)', {data(dd).latitude}, 'Un', 0));
    lon = cell2mat(cellfun(@(x) x(:)', {data(dd).longitude}, 'Un', 0)); lon = lon - 360;
    hs = cell2mat(cellfun(@(x) x(:)', {data(dd).hs}, 'Un', 0));
    swh_proc = cell2mat(cellfun(@(x) x(:)', {data(dd).swh_karin_proc}, 'Un', 0));
    swh_raw = cell2mat(cellfun(@(x) x(:)', {data(dd).swh_karin}, 'Un', 0));
    swh_qual = cell2mat(cellfun(@(x) x(:)', {data(dd).swh_karin_qual}, 'Un', 0));
    swh_nadir = cell2mat(cellfun(@(x) x(:)', {data(dd).swh_nadir_altimeter}, 'Un', 0));
    swh_model = cell2mat(cellfun(@(x) x(:)', {data(dd).swh_model}, 'Un', 0));
    swh_msk = cell2mat(cellfun(@(x1,x2) x1(:)' + x2(:)', {data(dd).swh_karin_proc}, {data(dd).swh_karin_uncorr_unsmoothed}, 'Un', 0));
    xtrk = cell2mat(cellfun(@(x) abs(x(:))', {data(dd).cross_track_distance}, 'Un', 0)); 
    proc = cellfun(@(x) x(:)', {data(dd).processing}, 'Un', 0); proc = [proc{:}];
    depth = cellfun(@(x,y) repmat(x, [1 size(y,2)]), {data(dd).depth},{data(dd).t}, 'Un', 0); depth = cell2mat(cellfun(@(x) x(:)', depth, 'Un', 0));
    disttocoast = cell2mat(cellfun(@(x) x(:)', {data(dd).distance_to_coast}, 'Un', 0));
    rainrate = cell2mat(cellfun(@(x) x(:)', {data(dd).rain_rate}, 'Un', 0));
    % allocate variables
    A = {x,y,z, t, lat, lon, hs, swh_nadir, swh_model, xtrk, proc, depth, disttocoast,swh_msk,swh_raw,swh_qual, swh_proc};

    % only keep data that overlaps in coverage
    nn = cellfun(@(fldnm) isnan(cell2mat(cellfun(@(x) x(:), {data(dd).(fldnm)}, 'Un', 0))), ogfldnms, 'Un', 0);
    nn = cell2mat(nn); nn = sum(nn,2); nn = find(~nn); 
    NN = nn; 
    A = cellfun(@(x) x(NN), A, 'Un', 0);
    [x,y,z,t, lat, lon, hs, swh_nadir, swh_model, xtrk, proc, depth, disttocoast,swh_msk,swh_raw,swh_qual, swh_proc] = A{:};



    % only processing as specified
    nn = ~cellfun(@(x) ~isempty(regexp(x,procname)), proc); nn = find(nn==0); 
    NN = nn;
    % only within good cross track region
    nn = find(isin(xtrk,[10 60])); 
    NN = intersect(NN,nn);
    % remove bad data (zero values)
    nn = find(swh_raw~=0); 
    NN = intersect(NN,nn);
    % remove outliars 
    nn = find(abs(swh_proc-hs)./hs < 1); 
    NN = intersect(NN,nn);


    % % final restrict then reset variable names
    A = cellfun(@(x) x(NN), A, 'Un', 0);
    [x,y,z, t, lat, lon, hs, swh_nadir, swh_model, xtrk, proc, depth, disttocoast,swh_msk,swh_raw,swh_qual, swh_proc] = A{:};


    % -------------------------------------------------------------------------
    
    if contains(fldnm,'swh_karin') & correctopt & ~contains(fldnm,'proc')
        corr = mdl(y, z); corr(isinf(corr)) = 0;
        nn = ~cellfun(@(x) ~isempty(regexp(x,'PIC2')), proc); nn = find(nn==0); 
        corr(nn) = 0; 
        corr(isinf(corr)) = NaN;
        y = y - corr;
    else
        corr = zeros(size(x));
    end

    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    setaxes(ha(hi),'YLabel', [ylb ' SWH [m]']);
    setaxes(ha(hi),'XLabel', [xlb ' $H_s$ [m]'] );
    clear sts
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    thisha = ha(hi); axes(thisha); hold on;

    scatter(x, y, 6, 'k', 'filled', 'MarkerFaceAlpha', 0.5)


    %%% ADD STATS LABELS
    dx = - 0.23; dy = 0.015;
    xdata = x; ydata = y; 
    sts.cc = ccfxn(xdata,ydata);
    sts.rmse = rmsefxn(xdata,ydata);
    sts.crmse = crmsefxn(xdata,ydata);
    sts.bias = biasfxn(xdata,ydata);
    x0 = sum(thisha.Position([1 3])) + dx; y0 = thisha.Position(2) + dy; 
    textbypos(x0,y0,...
        {   ['n$=$' num2str(sum(~isnan(xdata) & ~isnan(ydata)))], ...
            ['bias$=$' num2str(round(nanmean([sts.bias]),2)) 'm'], ...
            ['cRMSE$=$' num2str(round(nanmean([sts.crmse   ]),2)) 'm'],...
            ['CC$=$' num2str(round(nanmean([sts.cc]),2))]}, ...
        'interpreter', 'latex', 'VerticalAlignment','bottom');


    %%% PULL OUT HIGHLIGHT
    if hlopt
        tt = find(isin(t, thl_lims));
        if dd == 1;
            ll = find(isin(lat, [39 43]) & isin(lon, [-126 -123]));
            tt = intersect(tt,ll);
        end

        % scatter(xdata(tt), ydata(tt),12, 'r', 'filled') 
        scatter(x(tt), y(tt),51, x(tt)', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth',1, 'MarkerEdgeAlpha',0.5); 
        disp([xlb ' Hs range: ' num2str(range(x(tt)))])
        disp([ylb ' SWH range: ' num2str(range(y(tt)))]);
        
        % [xb yb] = drawbox(minmax(x(tt)), [0 10]);
        % patch(xb,yb,'k', 'FaceAlpha', 0.05, 'EdgeColor', 'none')
        % [xb yb] = drawbox([0 10],minmax(y(tt)));
        % patch(xb,yb,'k', 'FaceAlpha', 0.05, 'EdgeColor', 'none')
    end
     
    
end

quietaxes(ha(1),'x')
% AddLetters2Plots({ ha(1) ha(2)},...
%      {'(a)', '(b)', '(c)', '(d)', '(e)', '(f)', '(g)'},...
%     'BackgroundColor', 'w', 'Margin', 1,...
%     'HShift', 0.01, 'VShift', 0.015, ...
%     'FontName', 'Times', 'FontSize', 10, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top')





%% [save] do snapshot
dscl = 1.8;
dlat = 1.3; dlon = 1.26;
dlat = dlat*dscl; dlon = dlon*dscl;

dpopt = 1; 
sshopt = 0; 

if dd==2 
    % latlims = [34.8 36.6] + [0.9 -1.1].*0.25 ; lonlims = [-125.8 -124.2] + [1 -1].*0.17;
    latlims = 36 + [-1 1].*dlat./2;
    lonlims = -125.1 + [-1 1].*dlon./2;
    swot_fpath = [base_path '/data/SWOT/onedayrepeat/'];
else
    swot_fpath = [base_path '/data/SWOT/westcoast/'];
    
    % latlims = 41.3 + [-1 1].*(3.5/2.1); lonlims =  -124.9 + [-1 1].*(2.5/2.1);
    % latlims = latlims + [-1 1]*1;
    % lonlims = lonlims + [-1 1]*1;
    latlims = 41.1 + [-1 1].*dlat./2;
    lonlims = -124.9 + [-1 1].*dlon./2;
end
fpath = swot_fpath;

% -------------------------------------------------------------------------

fnames = dir([fpath '*.nc']);
ftimes = cellfun(@(x) x{8}, cellfun(@(x) strsplit(x, '_'), {fnames.name}, 'Un', 0), 'Un', 0);
ftimes = cellfun(@(x) datenum(x, 'yyyymmddTHHMMSS'), ftimes);
% tt = find(isin(ftimes, thl_lims)); fi = tt(1);
[dtss, ss] = sort(abs(ftimes-thl));
[dtf, fi] = min(abs(ftimes - thl));
swot = load_any_nc([fpath fnames(fi).name]);
% seconds since 2000-01-01 00:00:00.0
swot.time = swot.time./(3600*24) + datenum(2000,1,1);
% filter
Ny = size(swot.latitude,2);
yy = find(isin(swot.latitude(1,:) ,latlims + [-1 1].*1) & isin(swot.longitude(1,:) - 360, lonlims  + [-1 1].*1)); 
sfldnms = fieldnames(swot);
for ei = 1:length(sfldnms)
    fldnm = sfldnms{ei};
    tmp = swot.(fldnm);
    if size(tmp,1) == Ny
        swot.(fldnm) = tmp(yy,:);
    elseif size(tmp,2) == Ny
        swot.(fldnm) = tmp(:,yy); 
    end
end
swot.fname = fnames(fi).name; 
[swot.t0, swot.cycle, swot.pass, swot.processing, swot.orientation, swot.angle] = get_swot_info(swot, swot.fname);
modopt = 0; smoothopt = 1; interpopt = 1; patchopt = 0; Lavg = 5;  % anomaly patching on SWH
varargin = {'patch', patchopt, 'qcmodel', modopt, 'model', mdl, 'smooth', smoothopt, 'Lavg', Lavg, 'correct', 1, 'mask', 1, 'mindepth', -1, 'interp', interpopt};

        fldnm = 'swh_karin';
        Z = process_swot(swot, fldnm, varargin{:});
        swot.([fldnm '_proc']) = Z;
        Z = process_swot(swot, fldnm, varargin{:}, 'smooth', 0);
        swot.([fldnm '_proc_unsmoothed']) = Z;
    
        Z = process_swot(swot, fldnm, varargin{:}, 'patch', 2);
        swot.([fldnm '_proc_patch']) = Z;
        Z = process_swot(swot, fldnm, varargin{:}, 'smooth', 0, 'patch', 2);
        swot.([fldnm '_proc_patch_unsmoothed']) = Z;
    
        Z = process_swot(swot, fldnm, varargin{:}, 'smooth', 0,'patch',0, 'correct',0);
        swot.([fldnm '_uncorr_unsmoothed']) = Z;
        Z = process_swot(swot, fldnm, varargin{:}, 'smooth', 0,'patch',2, 'correct',0);
        swot.([fldnm '_uncorr_patch_unsmoothed']) = Z;
        Z = process_swot(swot, fldnm, varargin{:}, 'correct',0);
        swot.([fldnm '_uncorr']) = Z;     
        Z = process_swot(swot, fldnm, varargin{:},'correct',0,'patch',2);     
        swot.([fldnm '_uncorr_patch']) = Z;
        fldnm = 'ssha_karin_2';
        Z = process_swot(swot, fldnm, varargin{:}, 'patch', 2);
        swot.([fldnm '_proc_patch']) = Z;


% -------------------------------------------------------------------------
figure(405); clf; 
setfigsize(gcf, [352   627])
% tight_subplot(Nh, Nw, [gap_h gap_w], [lower upper], [left right]) 
ha = tight_subplot(2, 1, [0.02 0.02], [0.08 0.015], [0.13 0.03]);

% -------------------------------------------------------------------------
setaxes(ha(1:end),'XLim', lonlims);
setaxes(ha(1:end),'YLim', latlims);
setaxes(ha(:),'CLim', clims);
setaxes(ha(:),'Colormap', subsetcmap(nicejet, 100));
format_fig4print(ha)
setaxes(ha(:),'DataAspectRatio', [1 1 1]);
format_mapticks(ha, 'xy')
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


for hi=1:length(ha)
    disp(hi)
    dd = dds(hi);
    xlb = dlabs{hi};
    fldnm = fldnms{hi}; fldnm = strrep(fldnm, 'uncorr', 'proc');
    ylb = fldlbs{hi};

    t = cell2mat(cellfun(@(x) x(:)', {data(dd).t}, 'Un', 0));
    lat = cell2mat(cellfun(@(x) x(:)', {data(dd).latitude}, 'Un', 0));
    lon = cell2mat(cellfun(@(x) x(:)', {data(dd).longitude}, 'Un', 0)); lon = lon - 360;
     hs = cell2mat(cellfun(@(x) x(:)', {data(dd).hs}, 'Un', 0));
     dp = cell2mat(cellfun(@(x) x(:)', {data(dd).dp}, 'Un', 0));
    swh_proc = cell2mat(cellfun(@(x) x(:)', {data(dd).swh_karin_proc}, 'Un', 0));
    swh_raw = cell2mat(cellfun(@(x) x(:)', {data(dd).swh_karin}, 'Un', 0));
    xtrk = cell2mat(cellfun(@(x) abs(x(:))', {data(dd).cross_track_distance}, 'Un', 0)); 
    proc = cellfun(@(x) x(:)', {data(dd).processing}, 'Un', 0); proc = [proc{:}];
    

    A = {t,lon, lat, hs,  dp , swh_proc, swh_raw, xtrk,proc};

    % only keep data that overlaps in coverage
    nn = cellfun(@(fldnm) isnan(cell2mat(cellfun(@(x) x(:), {data(dd).(fldnm)}, 'Un', 0))), ogfldnms, 'Un', 0);
    nn = cell2mat(nn); nn = sum(nn,2); nn = find(~nn); 
    NN = nn; 
    A = cellfun(@(x) x(NN), A, 'Un', 0);
    [t,lon, lat, hs,  dp , swh_proc, swh_raw, xtrk,proc] = A{:};



    % only processing as specified
    nn = ~cellfun(@(x) ~isempty(regexp(x,procname)), proc); nn = find(nn==0); 
    NN = nn;
    % only within good cross track region
    nn = find(isin(xtrk,[10 60])); 
    NN = intersect(NN,nn);
    % remove bad data (zero values)
    nn = find(swh_raw~=0); 
    NN = intersect(NN,nn);
    % remove outliars 
    nn = find(abs(swh_proc-hs)./hs < 1); 
    NN = intersect(NN,nn);


    % % final restrict then reset variable names
    A = cellfun(@(x) x(NN), A, 'Un', 0);
    [t,lon, lat, hs,  dp , swh_proc, swh_raw, xtrk,proc] = A{:};



    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    X = swot.longitude - 360; Y = swot.latitude; 
    msk = ones(size(X)); msk(~isin(abs(swot.cross_track_distance./1000), [10 60])) = NaN;
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    thisha = ha(hi); axes(thisha); hold on;
    Z = swot.(fldnm);
    pcolor(X,Y,Z.*msk); shading flat;

    if sshopt
        Z = swot.ssha_karin_2_proc_patch + swot.height_cor_xover;
        contour(X,Y,Z.*msk,[-1:0.02:1], 'k-', 'EdgeAlpha',0.3)
    end
        

    if contains(fldnm, 'karin');
        xi = floor(size(X,1)/2);
        scatter(swot.longitude_nadir-360, swot.latitude_nadir, 30, swot.swh_nadir_altimeter(xi,:), 'filled')
    end

    if hlopt
        tt = find(isin(t, thl_lims));
        ll = find(isin(lat, latlims + [-1 1]) & isin(lon, lonlims + [-1 1]));
        tt = intersect(tt,ll);
        
        if dpopt
            scl = 0.13; 
            quiver(lon(tt), lat(tt), scl.*cosd(met2oc(oc2pol(dp(tt)))),scl.*sind(met2oc(oc2pol(dp(tt)))), 'k-', 'AutoScale', 'off' );
            disp(['MEAN WAVE DIRECTION FROM BUOYS'])
            disp(num2str(dirnanmean(dp(tt))))
        end
        scatter(lon(tt), lat(tt),80, hs(tt)', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth',1, 'MarkerEdgeAlpha',0.5); 



    end
    if dd~=2
        plot_mcoast(mcoast, lonlims, latlims, 'res',5);
    end


    c = fixedcolorbar(gca, 'Location','southoutside', 'FontSize',9, 'AxisLocationMode', 'manual', 'YAxisLocation', 'top'); 
    if dd==2; c.Ticks = [2:0.5:4]; end
    c.Position(3) = c.Position(3)/4;
    c.Position(4) = c.Position(4)*0.5;
    c.Position(1) = sum(thisha.Position([1 3])) - c.Position(3) - 0.075; 
    c.Position(2) = sum(thisha.Position([2])) + 0.015; 
    ylabel(c, [ylb ' SWH [m]'], 'Interpreter', 'latex', 'FontSize',8);
    % ylabel(c, ['SWH [m]'], 'Interpreter', 'latex', 'FontSize',8);
end



quietaxes(ha(1),'x')


