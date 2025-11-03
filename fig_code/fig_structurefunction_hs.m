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
    % cleanLegend.m
    % quietaxes.m
    % textbypos.m
    % setaxes.m
    % format_fig4print.m
    % subsetcmap.m

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


%% [save] one plot for moorings, structure function in units meters
fs = 15;

dd = 2;
bins = [1 5:10:100 125:25:200 300:100:1000]; 
nmin = 2;


hslims = [1 10];
clear sts
mipairs = nchoosek(1:length(data(dd).id),2);
distpairs = NaN(length(mipairs), size(data(dd).t,2));
for mi=1:length(mipairs) 
    mi1 = mipairs(mi,1); mi2 = mipairs(mi,2);
    % dist = distance(data(dd).lat0(mi1), data(dd).lon0(mi1), data(dd).lat0(mi2), data(dd).lon0(mi2), referenceEllipsoid('wgs84'))./1000;
    dist = distance(data(dd).lat(mi1,:), data(dd).lon(mi1,:), data(dd).lat(mi2,:), data(dd).lon(mi2,:), referenceEllipsoid('wgs84'))./1000;
    distpairs(mi,:) = dist;
end

 
qualflag = 1;
uncflag = 0;
corropt = 1;

% -------------------------------------------------------------------------
figure(69954); clf; 
setfigsize(gcf, [537         407])
% tight_subplot(Nh, Nw, [gap_h gap_w], [lower upper], [left right]) 
% oldha = tight_subplot(2, 4, [0.025 0.05], [0.11 0.05], [0.05 0.05]);
% [ha] = reorg_tightsubplot(oldha, {[1 2],[4 5]+1, [3 7], [4 8
ha = tight_subplot(1, 1, [0.025 0.015], [0.12 0.025], [0.12 0.03]);
% [ha] = reorg_tightsubplot(oldha, {[1 2],[4 5], [3 6]});

% -------------------------------------------------------------------------
setaxes(ha(:),'YLim', [0 1].*0.35);
setaxes(ha(:),'CLim', [0 70]);
% setaxes(ha(:),'Colormap', subsetcmap(nicejet, 20));
linkaxes(ha(:))

format_fig4print(ha, 'FontSize',fs)
setaxes(ha(:),'XLabel', 'separation distance $r$ [km]');
setaxes(ha(:),'YLabel', '$|SWH(x) - SWH(x+r)|$ [m]');
% quietaxes(ha(2:end), 'y');

setloglog(ha)
setaxes(ha(:),'XDir', 'reverse');
if dd == 1; xlims = [6 1100]; ylims = [0.000005 2]*0.1; else; xlims = [6 130]; ylims = [0.000005 0.23]*30; end;
ylims = [0.002 0.8]
setaxes(ha(:),'XLim', xlims);
setaxes(ha(:),'YLim', ylims);
setaxes(ha(:),'XTick', [10 20 50 100  200 500]);
setaxes(ha(:), 'YMinorGrid', 'off')

setaxes(ha(:),'CLim', [0 max(distpairs)]);
setaxes(ha(:),'CLim', [0.5 5]);

setaxes(ha(:),'YTick', [0.01 0.05 0.1 0.2]);

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
fldnms = { 'swh_karin_qc', 'swh_nadir_altimeter', 'swh_model'};
fldquas = { 'swh_karin_qual_avg', 'swh_karin_qual_avg', 'swh_karin_qual_avg'};
fldlbs = { 'SWOT KaRIn', 'SWOT nadir', 'model', 'SWOT KaRIn'};

fldnms = {  'swh_model', 'swh_karin_qc',};
fldquas = { 'swh_karin_qual_avg', 'swh_karin_qual_avg'};
fldlbs = {  'model', 'SWOT KaRIn',};


fldnms = {  'swh_model', 'swh_karin_proc_patch'}; 
% fldnms = {  'swh_model', 'swh_karin_proc'}; 
% fldnms = {  'swh_model', 'swh_karin_proc_unsmoothed'}; 
% fldnms = {  'swh_model','swh_karin_proc_patch', 'swh_karin_proc_unsmoothed'}; 
fldlbs = {  'model', 'SWOT KaRIn','SWOT KaRIn'};

cmap = subsetcmap(nicejet, length(fldnms)+4); cmap = flipud(cmap([end-4 end],:)); 
cmap(2,:) = cmap(2,:) + [-0.2 -0.15 0.15];
cmap = flipud(cmap);
cmap(cmap>1) = 1;
cmap = [cmap; cmap(2,:)]

for hi=1:length(fldnms)
    fldnm = fldnms{hi};
    ylb = fldlbs{hi};
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    t_swot = data(dd).t(1,:); 
    var_mo = data(dd).hs';
    var_swot = data(dd).(fldnm)';
    xtrck_swot = data(dd).cross_track_distance'; xtrck_swot = abs(xtrck_swot);
    NN = isnan(data(dd).swh_karin_uncorr' + data(dd).swh_model' + data(dd).swh_nadir_altimeter'); 
    NN = NN | squeeze(sum(cell2mat(cellfun(@(fldnm) permute(isnan(data(dd).(fldnm)), [1 3 2]), fldnms, 'Un', 0)),2))';
    var_mo(NN) = NaN; var_swot(NN) = NaN; 
    if corropt & contains(fldnm, 'karin') & ~contains(fldnm, 'proc')
       corr = mdl(var_swot, xtrck_swot);
       var_swot = var_swot - corr; 
    end

    xdata = NaN(length(t_swot), length(mipairs));
    ydata = NaN(length(t_swot), length(mipairs));
    zdata = NaN(length(t_swot), length(mipairs));


    for mi=1:length(mipairs)
        mi1 = mipairs(mi,1); mi2 = mipairs(mi,2);
        % dist = distpairs(mi);
        xdata(:,mi) = var_mo(:,mi1) - var_mo(:,mi2);
        ydata(:,mi) = var_swot(:,mi1) - var_swot(:,mi2);
        var0 = nanmean([var_mo(:,mi1) var_mo(:,mi2)]');
        zdata(:,mi) = var0;
    end

    bdata = distpairs';
    xdata = xdata(:); ydata = ydata(:); zdata = zdata(:); bdata = bdata(:);


    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    thisha = ha(1); axes(thisha); hold on; 
    ii = find(~isin(zdata, hslims)); xdata(ii) = NaN; ydata(ii) = NaN;
    nn = find(~isnan(xdata) & ~isnan(ydata)); xdata = xdata(nn); ydata = ydata(nn); zdata = zdata(nn); bdata = bdata(nn); 
    
    bindata = bdata; 
    n = 1; 
    if hi==1
        var = abs(xdata).^n;  var = abs(var); col = [0.8 0.1 0.2];  col = [0 0 0]+0.2;
        sts = boxplot_bybin(gca, var, bindata, bins, 'nmin',5, 'color', col, 'FaceAlpha',0.6, 'wid', 0.75.*mode(diff(bins)/2), 'plot', 0);
        nn = find(sts.n>=nmin);
        % shaded_limits(sts.x0(nn), sts.q2(nn), [sts.p10(nn)'; sts.p90(nn)' ], col, 'LineWidth',2, 'FaceAlpha', 0.05)
        shaded_limits(sts.x0(nn), sts.q2(nn), [sts.q1(nn)'; sts.q3(nn)' ], col, 'LineWidth',2, 'FaceAlpha', 0.1)
        scatter(sts.x0(nn), sts.q2(nn), 40, 'k', 'filled', 's', 'MarkerFaceColor', col, 'LineWidth',0.5, 'MarkerEdgeColor', 'k' , 'Displayname', 'GNSS buoy $H_s$')
        plot(sts.x0(nn), sts.q3(nn), 'k:', 'Color', col)
        plot(sts.x0(nn), sts.q1(nn), 'k:', 'Color', col)
    end

    var = abs(ydata).^n; var = abs(var); col = [0.2 0.1 0.8]; col = cmap(hi,:);
    sts = boxplot_bybin(gca, var, bindata, bins, 'nmin',5, 'color', col, 'FaceAlpha',0.6, 'wid', 0.75.*mode(diff(bins)/2), 'plot', 0);
    nn = find(sts.n>=nmin);
    % shaded_limits(sts.x0(nn), sts.q2(nn), [sts.p10(nn)'; sts.p90(nn)' ], col, 'LineWidth',2, 'FaceAlpha', 0.1)
    shaded_limits(sts.x0(nn), sts.q2(nn), [sts.q1(nn)'; sts.q3(nn)' ], col, 'LineWidth',2, 'FaceAlpha', 0.2)
    scatter(sts.x0(nn), sts.q2(nn), 40, 'k', 'filled', 's', 'MarkerFaceColor', col, 'LineWidth',0.5, 'MarkerEdgeColor', 'k' , 'Displayname', [ylb ' SWH'])
    plot(sts.x0(nn), sts.q3(nn), 'k:', 'Color', col)
    plot(sts.x0(nn), sts.q1(nn), 'k:', 'Color', col)

end

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
thisha = ha(1); axes(thisha); hold on; 
cleanLegend(gca, 'northeast', 'NumColumns',3, 'Interpreter', 'latex', 'FontSize', fs-3)



% savejpg(gcf, 'fig_structurefunction_hs', [base_path(1:12) 'desktop/'], 'on')