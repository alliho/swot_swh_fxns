base_path = '/Users/ajho/Documents/JPL/';
local_path = [base_path 'papers/swot_swh_calval/swot_swh_fxns/'];
addpath([local_path 'swh_fxns/matlab/'])
addpath([local_path 'fig_code/addfxns/'])
swot_fpath = [base_path '/data/SWOT/onedayrepeat/'];

%%% ADD DATA
load('/Users/ajho/Documents/myrepos/supportingdata/coastline_labeled.mat');

%%% OTHER DEPENDENCIES:
% m_map

%% LOAD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%d

%% load SWOT
clear SWOT
latlims = [32.3 43] + [-1 1]*180; lonlims = [-131.3 -118.2] + [-1 1]*180;

fpath = swot_fpath;
daterange = datenum(2023,5,1.5) + [-1 0];
SWOT = load_swot(fpath, lonlims, latlims, daterange);

fpath = '/Volumes/ahomiscdata/data_JPL/SWOT/misc/';
SWOT_back = load_swot(fpath, 'search', 'SWOT_L2_LR_SSH_Expert_577_014_20230709T170118_20230709T175142_PGC0_02.nc');
SWOT = [SWOT SWOT_back];
fpath = '/Volumes/ahomiscdata/data_JPL/SWOT/misc/';
SWOT_back = load_swot(fpath, 'search', 'SWOT_L2_LR_SSH_Expert_577_025_20230710T022320_20230710T031427_PGC0_02.nc');
SWOT = [SWOT SWOT_back];


%% load NDBC buoys
fpath = [local_path 'data/'];
load([fpath 'buoys_SWOT.mat']); 
buoys = data; 

%% load GNSS Cal/Val buoys
fpath = [local_path 'data/'];
fnames = dir([fpath 'GNSS*.nc']);
clear MO
for fi=1:length(fnames)
    mo = load_any_nc([fpath fnames(fi).name]);
    mo.t = mo.time./(3600*24) + mo.timeinfo.t0;
    MO(fi) = mo;
end

%% [save] summary map with 3 subsets and buoy image (global crossover and zoomed) AND BUOYS
latlims = [34.8 36.6] + [1 -1].*0.1 ; lonlims = [-125.8 -124.2] + [1 -1].*0.1;
t0 = datenum(2023,5,1.5); clims = [3.2 4.8];


[~, si] = min(abs(cellfun(@nanmean, {SWOT.time}) - t0));
swot = SWOT(si);

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

figure(1058); clf;  
setfigsize(gcf, [1168         312])
% tight_subplot(Nh, Nw, [gap_h gap_w], [lower upper], [left right]) 
ha = tight_subplot(1, 4, [0.015 0.04], [0.12 0.023], [0.07 0.03]);
% [ha] = reorg_tightsubplot(oldha, {[1 5], [2 6], [3:4], [7:8]});
% -------------------------------------------------------------------------
ha(1).Position(1) = ha(1).Position(1) - 0.035;
ha(1).Position(3) = ha(1).Position(3) + 0.025;
dx = 0.025;
% -------------------------------------------------------------------------
format_fig4print(ha)
setaxes(ha(2:4), 'DataAspectRatio', [1 1 1])
setaxes(ha(3), 'XLim', lonlims)
setaxes(ha(3), 'YLim', latlims)
setaxes(ha(2), 'XLim', lonlims + range(lonlims).*[-0.5 1.5]*3)
setaxes(ha(2), 'YLim', latlims + range(latlims).*[-1 1]*3)
format_mapticks(ha(2:3), 'x')
format_mapticks(ha(2:3), 'y')
swocol = [0.25 0.29 0.6];

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
thisha = ha(1); axes(thisha); hold on; 
m_proj('ortho','lat',nanmean(latlims)+4,'long',nanmean(lonlims)'+20);
m_coast('patch',[1 1 1].*0.7);
m_grid('linest','-','xticklabels',[],'yticklabels',[]);
m_scatter([buoys.lon0], [buoys.lat0], 7,  [1 0.85 0.3], 'filled', 'MarkerEdgeColor', 'k', 'MarkerEdgeAlpha', 0.5);
m_plot(cellfun(@nanmean, {MO.lon}),cellfun(@nanmean, {MO.lat}), 'k.', 'MarkerSize',4)

[xpts,ypts] = drawbox(lonlims + range(lonlims).*[-0.5 1.5]*3,latlims + range(latlims).*[-1 1]*3);
m_patch(xpts, ypts, 'w', 'EdgeColor', 'k', 'LineWidth',2, 'FaceColor', 'none');


for si=1:length(SWOT)
    swot = SWOT(si); 
    xi = 31;
    X = swot.longitude; Y = swot.latitude; X(X>360) = X(X>360)-360;
    X = X([5 xi],:); Y = Y([1 xi],:); X = [X(1,:) fliplr(X(end,:))]; Y = [Y(1,:) fliplr(Y(end,:))]; 
    m_patch(X,Y,swocol+0.1, 'FaceAlpha', 0.4, 'EdgeAlpha', 0.2)
    xi = 69-xi + 2;
    X = swot.longitude; Y = swot.latitude;  X(X>360) = X(X>360)-360;
    X = X([xi end-3],:); Y = Y([xi end],:); X = [X(1,:) fliplr(X(end,:))]; Y = [Y(1,:) fliplr(Y(end,:))]; 
    m_patch(X,Y,swocol+0.1, 'FaceAlpha', 0.4, 'EdgeAlpha', 0.2)
 end


for hi=2:3
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    thisha = ha(hi); axes(thisha); hold on; 
    
    plot_mcoast(mcoast, lonlims, latlims)
    
    for si=1:length(SWOT)
        swot = SWOT(si); 
        xi = 31;
        X = swot.longitude; Y = swot.latitude; X(X>180) = X(X>180)-360;
        X = X([5 xi],:); Y = Y([1 xi],:); X = [X(1,:) fliplr(X(end,:))]; Y = [Y(1,:) fliplr(Y(end,:))]; 
        patch(X,Y,swocol+0.1, 'FaceAlpha', 0.4, 'EdgeAlpha', 0.2)
        xi = 69-xi + 2;
        X = swot.longitude; Y = swot.latitude;  X(X>180) = X(X>180)-360;
        X = X([xi end-3],:); Y = Y([xi end],:); X = [X(1,:) fliplr(X(end,:))]; Y = [Y(1,:) fliplr(Y(end,:))]; 
        patch(X,Y,swocol+0.1, 'FaceAlpha', 0.4, 'EdgeAlpha', 0.2)
    end
    
    if hi==2
        plot([buoys.lon0], [buoys.lat0], 'ko', 'MarkerFaceColor', [1 0.85 0.3], 'MarkerSize',4);
        [xpts,ypts] = drawbox(lonlims,latlims);
        patch(xpts, ypts, 'w', 'EdgeColor', 'k', 'LineWidth',2, 'FaceColor', 'none');
    end

    for mi=1:length(MO)
        [~, ti] = min(abs(MO(mi).t  - nanmean(swot.time)));
        x = MO(mi).lon; y = MO(mi).lat; 
        % plot(x,y, 'k-');
        if hi==2
            plot(nanmean(x), nanmean(y), 'k.', 'LineWidth',2);
        else
            plot(x,y, 'k-');
        end
    end
end


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
thisha = ha(4); axes(thisha); hold on;
im = imread([local_path 'fig_code/pmel_gps_mooring.png']);
imagesc(im);
set(gca, 'YDir', 'reverse');
xlim([0 size(im,2)*0.92])
ylim([0 size(im,1)])
set(gca, 'DataAspectRatio', [1 1 1]);
% grid off; 
set(gca, 'XTick', [])
set(gca, 'YTick', [])
% -------------------------------------------------------------------------

AddLetters2Plots({ ha(1) ha(2) ha(3) ha(4)},...
     {'(a)', '(b)', '(c)', '(d)', '(e)', '(f)', '(g)'},...
    'BackgroundColor', 'none', 'Margin', 1,...
    'HShift', -0.01, 'VShift', -0.035, ...
    'FontName', 'Times', 'FontSize', 10, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', 'Location', 'SouthEast')

