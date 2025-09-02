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
% -------------------------------------------------------------------------
%%% compute model
mdl = correct_swotswh('data', 'moorings');
mdl = correct_swotswh('data', 'buoys');
mdl = correct_swotswh();

%% [save] BIAS performance vs hs and cross track all moorings together | with histogram separaterd

addmdl = 0;

fldnm = 'swh_karin_qc_avg'
fldnm = 'swh_karin_qc'
fldnm = 'swh_karin'
fldnm = 'swh_karin_uncorr_patch';
% fldnm = 'swh_karin_uncorr_patch_unsmoothed'
% fldnm = 'swh_nadir_altimeter'
% fldnm = 'swh_model'
dd = 2; nmin = 5;
% dd = 1:2; nmin = 10;
dd = [1:3]; nmin = 15;
% dd = [1 3]; nmin = 15;

procname = 'P..0'; nmin = 20;
histlims = [0 800]; histtx = [250 500 750];
% procname = 'PIC2'; nmin = 20;
% histlims = [0 160]; histtx = [50 100 150];


hs = cell2mat(cellfun(@(x) x(:)', {data(dd).hs}, 'Un', 0));
swh = cell2mat(cellfun(@(x) x(:)', {data(dd).(fldnm)}, 'Un', 0));
swh_nadir = cell2mat(cellfun(@(x) x(:)', {data(dd).swh_nadir_altimeter}, 'Un', 0));
swh_model = cell2mat(cellfun(@(x) x(:)', {data(dd).swh_model}, 'Un', 0));
swh_msk = cell2mat(cellfun(@(x1,x2) x1(:)' + x2(:)', {data(dd).swh_karin_proc}, {data(dd).swh_karin_uncorr}, 'Un', 0));
xtrk = cell2mat(cellfun(@(x) abs(x(:))', {data(dd).cross_track_distance}, 'Un', 0)); 
proc = cellfun(@(x) x(:)', {data(dd).processing}, 'Un', 0); proc = [proc{:}];
swh_corr = mdl(swh, xtrk);
% 
% 
nn = find(~cellfun(@isempty, proc));
A = {hs, swh, swh_nadir, swh_model, xtrk, proc, swh_corr,swh_msk};
A = cellfun(@(x) x(nn), A, 'Un', 0);
[hs, swh, swh_nadir, swh_model, xtrk, proc, swh_corr,swh_msk] = A{:};




var = swh - hs; 
nn = find(cellfun(@(x) strcmp(x, 'PIC0') | strcmp(x, 'PGC0'), proc)); 
var(nn) = var(nn) - swh_corr(nn); 
var = var./hs;
var = abs(var);
q1 = prctile(var, [25]); q2 = prctile(var, [50]); q3 = prctile(var, [75]); iqr = q3-q1;
lims = q2 + 2.*[-1 1].*iqr;
NN = isnan(swh) | isnan(hs) | isnan(swh_nadir) | isnan(swh_model) | isnan(xtrk) | isnan(swh_msk);
NN = NN + ~isin(var, lims); 
NN = NN + ~cellfun(@(x) ~isempty(regexp(x,procname)), proc);
NN = NN ==0; 

NN = find(NN); 
A = {hs, swh, swh_nadir, swh_model, xtrk, proc, swh_corr,swh_msk};
A = cellfun(@(x) x(NN), A, 'Un', 0);
[hs, swh, swh_nadir, swh_model, xtrk, proc, swh_corr,swh_msk] = A{:};

xdata = hs; ydata = xtrk; 
xdata = swh;
zdata = swh - hs; %zdata = zdata./hs; 
pdata = proc; 

% -------------------------------------------------------------------------

cmap = subsetcmap(nicejet, 6); cmap = cmap(2:2:end,:);
cmap = [nicejet(1,:); 0 0 0; nicejet(end,:)];
cmap = [0.65 0.15 0.25; 0.3 0.3 0.3; 0.35 0.55 0.2]+0.1;
size(cmap)
% -------------------------------------------------------------------------

figure(2016); clf; 
setfigsize(gcf, [650         255])
% tight_subplot(Nh, Nw, [gap_h gap_w], [lower upper], [left right]) 
ha = tight_subplot(2, 2, [0.0 0.023], [0.18 0.025], [0.089 0.025]);
% setfigsize(gcf, [775         275])
% % tight_subplot(Nh, Nw, [gap_h gap_w], [lower upper], [left right]) 
% ha = tight_subplot(1, 2, [0.015 0.023], [0.16 0.025], [0.08 0.025]);
% [ha] = reorg_tightsubplot(oldha, {[1],[2:3]});
% -------------------------------------------------------------------------
for hi=[1 2]
    ha(hi).Position(2) = ha(hi+2).Position(2); 
    ha(hi).Position(4) =  ha(hi).Position(4)*2; 
end
for hi=[3 4]
    ha(hi).Position(4) =  ha(hi).Position(4)/4; 
end
% -------------------------------------------------------------------------
format_fig4print(ha)
setaxes(ha, 'FontSize', 16)
setaxes(ha(1:2), 'YLim', [-1 0.4])
setaxes(ha(1:2), 'YLim', [-1 1].*1 - 0.2)
setaxes(ha([1 3]), 'XLim', [0 6.35])
setaxes(ha([1 3]), 'XLim', [0 5.35])
setaxes(ha([1 3]), 'XLim', [-0.2 4.7])
setaxes(ha([1 3]), 'XLabel', 'SWH [m]' )
setaxes(ha(1:2), 'YLabel', '(KaRIn SWH - $H_s$) [m]' )
setaxes(ha([2 4]), 'XLim', [3 68])
setaxes(ha([2 4]), 'XLabel', 'x-track distance [km]')
prcs = [10 90];
prcs = [1 99];
facecol = [0.28 0.33 .5]+0.3
linkaxes(ha(1:2), 'y')


% -------------------------------------------------------------------------
bins = [-1:0.5:7.5]; 
% bins = [-0.25:0.5:5.25]; 
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
thisha = ha(1); axes(thisha); hold on; 
% errorbar(bins(nms)+ mode(diff(bins))/2, binmean, cellfun(@(x) x(1), binprcs), cellfun(@(x) x(2), binprcs), 'k-',...
%     'LineWidth', 2, 'LineStyle','-', 'Color', 'k');
% plot(bins(nms)+ mode(diff(bins))/2, binmean, 'k.', 'LineWidth',2, 'MarkerSize',20, 'Color', 'k', 'DisplayName', num2str(mi));
% yline(0)
% text(bins(nms) + mode(diff(bins))/2, 0.3*ones(size(bins(nms))),num2str(binN'), 'HorizontalAlignment', 'center', 'FontName', 'times', 'Color', [1 1 1].*0.8)

[sts] = boxplot_bybin(gca, zdata, xdata, bins, 'nmin',nmin, 'color', facecol, 'FaceAlpha',1, 'wid', 0.75.*mode(diff(bins)/2));
nn = find(sts.n>=nmin);
scatter(sts.x0(nn), sts.q2(nn), 40, 'k', 'filled', 's', 'MarkerFaceColor', facecol, 'LineWidth',1, 'MarkerEdgeColor', 'k' , 'Displayname', 'all')
% text(sts.x0(nn), ones(size(nn)).*0.34, num2str(sts.n(nn)), 'HorizontalAlignment','center')
% scatter(sts.x0, sts.q2, 40, 'k', 'filled', 's', 'MarkerFaceColor', 'w', 'LineWidth',1, 'MarkerEdgeColor', 'k' )

% lw = 1
% errorbar(sts.x0,sts.q2,range([sts.q2; sts.q1-1.5.*sts.iqr]),range( [sts.q2; sts.q3+1.5.*sts.iqr]), 'r-', 'LineStyle', 'none', 'LineWidth',lw);
% errorbar(sts.x0,sts.q2,range([sts.q2; sts.q1]),range( [sts.q2; sts.q3]), 'r-', 'LineStyle', 'none', 'LineWidth',lw);
% errorbar(sts.x0,sts.q2,sts.q2 - sts.q1,sts.q3 - sts.q2, 'r-', 'LineStyle', 'none', 'LineWidth',lw);
% errorbar(sts.x0,sts.q2,sts.q2 - (sts.q1 - 1.5.*sts.iqr),sts.q2 - (sts.q3 - 1.5.*sts.iqr), 'r-', 'LineStyle', 'none', 'LineWidth',lw);
% errorbar(X,Y,YNEG,YPOS,XNEG,XPOS)


ii = find(ydata<30); col = cmap(end,:);
[sts] = boxplot_bybin(gca, zdata(ii), xdata(ii), bins, 'nmin',nmin, 'color', [0.4 0.6 0.3]+0.2, 'FaceAlpha',1, 'wid', 0.75.*mode(diff(bins)/2), 'plotopt',0);
nn = find(sts.n>=nmin);
plot(sts.x0(nn), sts.q2(nn), 'k:', 'Color', col)
scatter(sts.x0(nn), sts.q2(nn), 40, 'k', 'filled', '<', 'MarkerFaceColor', col, 'LineWidth',1, 'MarkerEdgeColor', 'k' , 'DisplayName', 'x-track $ < 30 $km')
% errorbar(sts.x0(nn), sts.q2(nn),range([sts.q2(nn); sts.q2(nn)-1.5.*sts.iqr(nn)]),range([sts.q2(nn); sts.q2(nn)+1.5.*sts.iqr(nn)]), 'k-', 'LineStyle', 'none', 'LineWidth',1, 'Color', col);
% shaded_limits(sts.x0(nn), sts.q2(nn), [sts.q2(nn)-1.5.*sts.iqr(nn);  sts.q2(nn)+1.5.*sts.iqr(nn)], col)
% shaded_limits(sts.x0(nn), sts.q2(nn), [sts.q1(nn);  sts.q3(nn)], col, 'FaceAlpha', 0.15, 'LineWidth',2, 'LineStyle', ':');

ii = find(ydata>=30); col = cmap(1,:);
[sts] = boxplot_bybin(gca, zdata(ii), xdata(ii), bins, 'nmin',nmin, 'color', col, 'FaceAlpha',2, 'wid', 0.75.*mode(diff(bins)/2), 'plotopt',0);
nn = find(sts.n>=nmin);
plot(sts.x0(nn), sts.q2(nn), 'k:', 'Color', col)
scatter(sts.x0(nn), sts.q2(nn), 40, 'k', 'filled', '>', 'MarkerFaceColor', col, 'LineWidth',1, 'MarkerEdgeColor', 'k' , 'DisplayName', 'x-track $ > 30 $km')
% errorbar(sts.x0(nn), sts.q2(nn),range([sts.q2(nn); sts.q2(nn)-1.5.*sts.iqr(nn)]),range([sts.q2(nn); sts.q2(nn)+1.5.*sts.iqr(nn)]), 'k-', 'LineStyle', 'none', 'LineWidth',1, 'Color', col);
% shaded_limits(sts.x0(nn), sts.q2(nn), [sts.q2(nn)-1.5.*sts.iqr(nn);  sts.q2(nn)+1.5.*sts.iqr(nn)], col)
% shaded_limits(sts.x0(nn), sts.q2(nn), [sts.q1(nn);  sts.q3(nn)], col, 'FaceAlpha', 0.15, 'LineWidth',2, 'LineStyle', ':');

% col = cmap(2,:);
% [sts] = boxplot_bybin(gca, zdata, xdata, bins, 'nmin',5, 'color', facecol, 'FaceAlpha',1, 'wid', 0.75.*mode(diff(bins)/2), 'plotopt',0);
% nn = find(sts.n>3);
% shaded_limits(sts.x0(nn), sts.q2(nn), [sts.q1(nn);  sts.q3(nn)], col, 'FaceAlpha', 0.15, 'LineWidth',2, 'LineStyle', '-');


xticks([1:5])

yline(0, 'k--')
leg = cleanLegend(gca, 'northeast', 'Interpreter', 'latex', 'FontSize',10, 'NumColumns', 3);
leg.ItemTokenSize(1) = leg.ItemTokenSize(1)/2;
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
thisha = ha(3); axes(thisha); hold on; 
histogram(xdata, bins, 'FaceColor', facecol*0.9, 'EdgeColor', 'none', 'Normalization','count', 'FaceAlpha',0.3);
set(gca, 'Color', 'none', 'XColor', 'none', 'YColor', 'none', 'XTick', [], 'YTick', [], 'YLim', histlims)
tx = histtx; 
yline(tx, 'k:', 'Color', [facecol*0.7 0.3])
text(ones(size(tx))*4.58, tx, num2str(tx'), 'HorizontalAlignment','right', 'VerticalAlignment','middle', 'FontName', 'Times', 'FontSize',8, 'Color', [facecol*0.8 0.3])

% -------------------------------------------------------------------------
bins = [0:5:70]; 
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
thisha = ha(2); axes(thisha); hold on; 

[sts] = boxplot_bybin(gca, zdata, ydata, bins, 'nmin',nmin, 'color', facecol, 'FaceAlpha',1, 'wid', 0.75.*mode(diff(bins)/2));
nn = find(sts.n>=nmin);
scatter(sts.x0(nn), sts.q2(nn), 40, 'k', 'filled', 's', 'MarkerFaceColor', facecol, 'LineWidth',1, 'MarkerEdgeColor', 'k' , 'Displayname', 'all')
% text(sts.x0(nn), ones(size(nn)).*0.34, num2str(sts.n(nn)), 'HorizontalAlignment','center')

ii = find(xdata>2);  col = [0.3 0.5 0.2]+0.2;
[sts] = boxplot_bybin(gca, zdata(ii), ydata(ii), bins, 'nmin',nmin, 'color', [0.4 0.6 0.3]+0.2, 'FaceAlpha',1, 'wid', 0.75.*mode(diff(bins)/2), 'plotopt',0);
nn = find(sts.n>=nmin);
% scatter(sts.x0(nn), sts.q2(nn), 40, 'k', 'filled', '^', 'MarkerFaceColor', 'w', 'LineWidth',1, 'MarkerEdgeColor', 'k' )
% shaded_limits(sts.x0(nn), sts.q2(nn), [sts.q1(nn);  sts.q3(nn)], col, 'FaceAlpha', 0.1, 'LineWidth',1, 'LineStyle', ':');
plot(sts.x0(nn), sts.q2(nn), 'k:', 'Color', col)
scatter(sts.x0(nn), sts.q2(nn), 40, 'k', 'filled', '^', 'MarkerFaceColor', col, 'LineWidth',1, 'MarkerEdgeColor', 'k' , 'DisplayName', 'SWH $> 2$m')

ii = find(xdata<2);col = [0.6 0.2 0.2]+0.2;
[sts] = boxplot_bybin(gca, zdata(ii), ydata(ii), bins, 'nmin',nmin, 'color', [0.4 0.6 0.3]+0.2, 'FaceAlpha',1, 'wid', 0.75.*mode(diff(bins)/2), 'plotopt',0);
nn = find(sts.n>=nmin);
% scatter(sts.x0(nn), sts.q2(nn), 40, 'k', 'filled', 'v', 'MarkerFaceColor', 'w', 'LineWidth',1, 'MarkerEdgeColor', 'k' )
% shaded_limits(sts.x0(nn), sts.q2(nn), [sts.q1(nn);  sts.q3(nn)], col, 'FaceAlpha', 0.1, 'LineWidth',1, 'LineStyle', ':');
plot(sts.x0(nn), sts.q2(nn), 'k:', 'Color', col)
scatter(sts.x0(nn), sts.q2(nn), 40, 'k', 'filled', 'v', 'MarkerFaceColor', col, 'LineWidth',1, 'MarkerEdgeColor', 'k' , 'DisplayName', 'SWH $ < 2$m')




% xticks([1:5])
yline(0, 'k--')


leg = cleanLegend(gca, 'northeast', 'Interpreter', 'latex', 'FontSize',10, 'NumColumns', 3);
leg.ItemTokenSize(1) = leg.ItemTokenSize(1)/2;
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
thisha = ha(4); axes(thisha); hold on; 
histogram(ydata, bins, 'FaceColor', facecol*0.9, 'EdgeColor', 'none', 'Normalization','count', 'FaceAlpha',0.3);
set(gca, 'Color', 'none', 'XColor', 'none', 'YColor', 'none', 'XTick', [], 'YTick', [], 'YLim', histlims)
tx = histtx; 
yline(tx, 'k:', 'Color', [facecol*0.7 0.3])
text(ones(size(tx))*66, tx, num2str(tx'), 'HorizontalAlignment','right', 'VerticalAlignment','middle', 'FontName', 'Times', 'FontSize',8, 'Color', [facecol*0.8 0.3])

if addmdl
    axes(ha(1));
    bins = [0.5:0.25:4.5]; 
    xg = bins; 
    
    
    ii = find(ydata<30); col = cmap(end,:);
    yg = mdl(xg,nanmedian(ydata(ii)));
    plot(xg,yg,'k--', 'Color', col, 'LineWidth',1);

    ii = find(ydata>=30);  col = cmap(1,:);
    yg = mdl(xg,nanmedian(ydata(ii)));
    plot(xg,yg,'k--', 'Color', col, 'LineWidth',1);


    axes(ha(2));
    bins = [10:5:60]; 
    xg = bins; 
    
    
    ii = find(xdata>2);  col = [0.3 0.5 0.2]+0.2;
    yg = mdl(nanmedian(xdata(ii)),xg);
    plot(xg,yg,'k--', 'Color', col, 'LineWidth',1);

    ii = find(xdata<2);col = [0.6 0.2 0.2]+0.2;
    yg = mdl(nanmedian(xdata(ii)),xg);
    plot(xg,yg,'k--', 'Color', col, 'LineWidth',1);

end



% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
quietaxes(ha(2), 'y')

AddLetters2Plots({ ha(1) ha(2)},...
     {'(a)', '(b)', '(c)'},...
    'BackgroundColor', 'w', 'Margin', 1,...
    'HShift', 0.01, 'VShift', 0.03, ...
    'FontName', 'Times', 'FontSize', 10, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top')

