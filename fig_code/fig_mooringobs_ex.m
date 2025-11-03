% when located in working directory (/swot_swh_fxns/fig_code/
run set_env.m

%%% OTHER DEPENDENCIES:
% m_map
% gitsio/code_universal/*/
    % isin.m
    % plot_mcoast.m
    % load_any_nc.m
    % setfigsize.m
    % drawbox.m
% tight_subplot.m
% AddLetters2Plots.m

%% LOAD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%d
%% load SWOT
clear SWOT
latlims = [32.3 43]; lonlims = [-131.3 -118.2];

fpath = swot_fpath;
daterange = [datenum(2023,3,14) datenum(2023,7,15)]
SWOT = load_swot(fpath, lonlims, latlims, daterange);

%% process SWOT

modopt = 0; smoothopt = 1; interpopt = 1; patchopt = 0; Lavg = 5;  % anomaly patching on SWH
varargin = {'patch', patchopt, 'smooth', smoothopt, 'Lavg', Lavg, 'mask', 1, 'mindepth', -1, 'interp', interpopt};

for si=1:length(SWOT)
    swot = SWOT(si);
    disp(si)

    tic

    fldnm = 'swh_karin';

    Z = process_swot(swot, fldnm, varargin{:}, 'smooth', 0,'patch',0, 'correct',0);
    SWOT(si).([fldnm '_uncorr_unsmoothed']) = Z;
    Z = process_swot(swot, fldnm, varargin{:}, 'smooth', 0,'patch',2, 'correct',0);
    SWOT(si).([fldnm '_uncorr_patch_unsmoothed']) = Z;
    Z = process_swot(swot, fldnm, varargin{:}, 'correct',0);
    SWOT(si).([fldnm '_uncorr']) = Z;     
    Z = process_swot(swot, fldnm, varargin{:},'correct',0,'patch',2);     
    SWOT(si).([fldnm '_uncorr_patch']) = Z;

    toc
end
disp('SWOT has been processed')


%% load GNSS Cal/Val buoys
fpath = [local_path 'data/'];
fnames = dir([fpath 'GNSS*.nc']);
clear MO
for fi=1:length(fnames)
    mo = load_any_nc([fpath fnames(fi).name]);
    mo.t = mo.time./(3600*24) + mo.timeinfo.t0;
    MO(fi) = mo;
end

%% co-locate GPS moorings to KaRIn 
interpopt = 1; % interpolate in 3x3 grid to buoy location
for mi=1:length(MO)
    disp(mi); mo = MO(mi); if allnan(mo.hs); continue; end
    tic
    mosw = colocate_swot(SWOT, mo.t, mo.lon, mo.lat, 'interp', interpopt);
    toc
    fldnms = fieldnames(mosw);
    for fi=1:length(fldnms)
        fldnm = fldnms{fi};
        MO(mi).swot.(fldnm) = mosw.(fldnm);
    end
end

%% [save] plot
fs = 15;
mi = 3; 
swotopt = 1; nadiropt = 0;
fldnm = 'swh_karin_uncorr_patch';

mo = MO(mi);
figure(306); clf;
setfigsize(gcf, [1293         304])
% tight_subplot(Nh, Nw, [gap_h gap_w], [lower upper], [left right]) 
ha = tight_subplot(2, 1, [0.025 0.05], [0.09 0.05], [0.06 0.07]);
% [ha] = reorg_tightsubplot(oldha, {[1],[2:3]});
% -------------------------------------------------------------------------
dy = 0.08;
ha(1).Position(2) = ha(1).Position(2) + dy; 
ha(1).Position(4) = ha(1).Position(4) - dy; 
ha(2).Position(4) = ha(2).Position(4) + dy; 

% -------------------------------------------------------------------------
format_fig4print(ha, 'FontSize',fs);
setaxes(ha, 'XLim', minmax(mo.t));
setaxes(ha, 'XLim', minmax(mo.swot.time) + [-1 1].*15);
% setaxes(ha, 'XLim', minmax(mo.swot.time) + [-1 1].*2);
setaxes(ha, 'XLim', minmax(mo.swot.time) + [-1 0].*2);

cmap = [ buildcmap([ nicejet(1,:).*0.7;nicejet(1,:)],15); nicejet; buildcmap([nicejet(end,:); nicejet(end,:).*0.7],23) ];
setaxes(ha, 'Colormap', cmap);
% colormap(gca,[nicejet; buildcmap([nicejet(end,:); nicejet(end,:).*0.7],23) ]);
setloglog(ha(1), 'y')
setaxes(ha(1), 'YLim', [0.045 0.5])
setaxes(ha(1), 'YTick', [0.05 0.1 0.2 0.5])
setaxes(ha(1), 'CLim', [-2.6 1.5])
setaxes(ha(1), 'YLabel', '$f$ [Hz]')
setaxes(ha(2), 'YLabel', '$H_s$ [m]')
setaxes(ha(2), 'YLim', [0 6.7])

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
thisha = ha(1); axes(thisha); hold on; 
nn = find(~isnan(mo.hs));
pcolor(mo.t(nn), mo.f, log10(mo.sf(:,nn))); shading flat;
c = fixedcolorbar(gca); ylabel(c, '$s(f)$ [m$^2$ s]', 'interpreter', 'latex')
autodatetick(gca, 'x')
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
thisha = ha(2); axes(thisha); hold on; 
nn = find(~isnan(mo.hs));
plot(mo.t(nn), mo.hs(nn), 'k-', 'DisplayName', 'GNSS $H_s$', 'LineWidth',2)
autodatetick(gca, 'x')
if swotopt;
    if nadiropt
        plot(mo.swot.time, mo.swot.swh_nadir_altimeter, 'c.', 'MarkerSize',9, 'DisplayName' ,'SWOT nadir SWH', 'Color',  [0.11 0.28 0.8]+[0.3 0.3 0.2]);
    end
    nn = find(~isnan(mo.swot.(fldnm))); 
    plot(mo.swot.time(nn), mo.swot.(fldnm)(nn), 'r.', 'MarkerSize',10, 'DisplayName' ,'SWOT KaRIn SWH')
end
cleanLegend(gca, 'northeast', 'NumColumns', 3, 'interpreter', 'latex', 'Box', 'off', 'FontSize', fs-3)
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
xlims = get(gca, 'XLim'); 
alignyaxes(ha, xlims(1)-4.5)
AddLetters2Plots({ ha(1) ha(2)},...
   {'(a)', '(b)', '(c)', '(d)', '(e)', '(f)', '(g)'},...
  'BackgroundColor', 'none', 'Margin', 1,...
  'HShift', 0.005, 'VShift', 0.015, ...
  'FontName', 'Times', 'FontSize', fs, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top')

quietaxes(ha(1), 'x')

% savejpg(gcf, 'fig_mooringobs_ex', [base_path(1:12) 'desktop/'], 'on')