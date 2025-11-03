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



%% [save for paper] just HS and DIST
fs = 15; 
[mdl] = correct_swotswh('data', 'moorings');

fldnm = 'swh_karin_qc_avg'
fldnm = 'swh_karin_qc'
fldnm = 'swh_karin'
fldnm = 'swh_karin_proc'
% fldnm = 'swh_karin_uncorr_patch_unsmoothed'
fldnm = 'swh_karin_uncorr_unsmoothed'
% fldnm = 'swh_karin_uncorr_unsmoothed'
% fldnm = 'swh_nadir_altimeter'
% fldnm = 'swh_model'
dd = 1; nmin = 5;
% dd = 2; nmin = 10; data(2).depth = isnan(data(2).depth) + 6000;
% dd = [1 3]; nmin = 2;

facecol = [0.28 0.33 .5]+0.3
procname = 'P..0';  facecol = [0.76 0.02 0.02];
procname = 'P..0';  facecol = [0.76 0.02 0.02];
procname = 'PIC2'; facecol = [0.1 0.5 1]-0.1;


hs = cell2mat(cellfun(@(x) x(:)', {data(dd).hs}, 'Un', 0));
dp = cell2mat(cellfun(@(x) x(:)', {data(dd).dp}, 'Un', 0));
tp = cell2mat(cellfun(@(x) x(:)', {data(dd).tp}, 'Un', 0));
swh_raw = cell2mat(cellfun(@(x) x(:)', {data(dd).swh_karin}, 'Un', 0));
swh_qual = cell2mat(cellfun(@(x) x(:)', {data(dd).swh_karin_qual}, 'Un', 0));
swh_unc = cell2mat(cellfun(@(x) x(:)', {data(dd).swh_karin_uncert}, 'Un', 0));
swh = cell2mat(cellfun(@(x) x(:)', {data(dd).(fldnm)}, 'Un', 0));
swh_nadir = cell2mat(cellfun(@(x) x(:)', {data(dd).swh_nadir_altimeter}, 'Un', 0));
swh_model = cell2mat(cellfun(@(x) x(:)', {data(dd).swh_model}, 'Un', 0));
xtrk = cell2mat(cellfun(@(x) abs(x(:))', {data(dd).cross_track_distance}, 'Un', 0)); 
proc = cellfun(@(x) x(:)', {data(dd).processing}, 'Un', 0); proc = [proc{:}];
depth = cellfun(@(x,y) repmat(x, [1 size(y,2)]), {data(dd).depth},{data(dd).t}, 'Un', 0); depth = cell2mat(cellfun(@(x) x(:)', depth, 'Un', 0));
disttocoast = cell2mat(cellfun(@(x) x(:)', {data(dd).distance_to_coast}, 'Un', 0));
rainflag = cell2mat(cellfun(@(x) x(:)', {data(dd).rain_flag}, 'Un', 0));
rainrate = cell2mat(cellfun(@(x) x(:)', {data(dd).rain_rate}, 'Un', 0));

swh_corr = mdl(swh, xtrk);
% 
% 
% nn = find(~cellfun(@isempty, proc) & ~isnan(hs));
nn = find(~cellfun(@isempty, proc) & ~isnan(hs) & ~isnan(swh_raw));
A = {hs,tp,dp, swh, swh_unc, swh_nadir, swh_model, xtrk, proc, swh_corr, depth, disttocoast,swh_raw,swh_qual,rainflag,rainrate};
A = cellfun(@(x) x(nn), A, 'Un', 0);
[hs,tp,dp, swh, swh_unc, swh_nadir, swh_model, xtrk, proc, swh_corr, depth, disttocoast,swh_raw,swh_qual,rainflag,rainrate] = A{:};


NN = find(isin(abs(xtrk), [10 60]));
A = {hs,tp,dp,swh, swh_unc, swh_nadir, swh_model, xtrk, proc, swh_corr, depth, disttocoast,swh_raw,swh_qual,rainflag,rainrate};
A = cellfun(@(x) x(NN), A, 'Un', 0);
[hs,tp,dp, swh, swh_unc, swh_nadir, swh_model, xtrk, proc, swh_corr, depth, disttocoast,swh_raw,swh_qual,rainflag,rainrate] = A{:};

disp(['percent usable: ' num2str(100*sum(~isnan(swh))./sum(~isnan(hs)))])

figure(6234); clf; 
setfigsize(gcf, [897         501])
% tight_subplot(Nh, Nw, [gap_h gap_w], [lower upper], [left right]) 
ha = tight_subplot(2, 2, [0.04 0.018], [0.1 0.035], [0.1 0.03]);
% [ha] = reorg_tightsubplot(oldha, {[1],[2:3]});
dy = 0.06; 
for hi=1:2
    ha(hi).Position([2 4]) = ha(hi).Position([2 4]) + [1 -1]*dy;
end
for hi=3:4
    ha(hi).Position([ 4]) = ha(hi).Position([ 4]) + [1]*dy;
end
% -------------------------------------------------------------------------
format_fig4print(ha, 'FontSize',fs);
setaxes(ha([1:2]), 'YLim', [0 1100])
setaxes(ha([1:2]), 'YLabel', 'N')
setaxes(ha([3:4]), 'YLim', [0 100])
setaxes(ha([3:4]), 'YLabel', '\%')
setaxes(ha(3:4), 'YTick', [0:25:100])
setloglog(ha([2 4]-1), 'x')
setaxes(ha([2 4]-1), 'XLim', [0.5 900])
setaxes(ha([2 4]-1), 'XLabel', 'distance from coast [km]')
setaxes(ha([2 4]-1), 'XTick', [2 5 10 25 100 500 1000])

setaxes(ha([1 3]+1), 'XLim', [0 5.5])
setaxes(ha([1 3]+1), 'XLabel', '$H_s$ [m]')
% setaxes(ha([2 4]), 'XTick', [0.5 5 10 25 100 500 1000])

quietaxes(ha(setdiff(1:length(ha), [length(ha)-1:length(ha)])), 'x')
quietaxes(ha(setdiff(1:length(ha), [1:2:end])), 'y')
% 
% % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
thisha = ha(2); axes(thisha); hold on; 
bins = [0:0.3:6]; binvar = hs;
% yyaxis left; set(gca, 'YColor', 'k'); 
histogram(binvar(~isnan(hs)),bins, 'FaceColor', 'k', 'DisplayName', 'coastal buoy', 'EdgeColor', 'none')
y0 = histcounts(binvar(~isnan(hs)),bins);
histogram(binvar(~isnan(swh) & ~isnan(hs)), bins, 'FaceColor',  [0.3 0.75 0.9], 'DisplayName', 'SWOT KaRIn', 'EdgeColor', 'none')
y = histcounts(binvar(~isnan(swh) & ~isnan(hs)),bins);
% yyaxis right; col = [0.7 0.1 0.2]; set(gca, 'YColor', col); ylim([0 100])
% x = bins(1:end-1) + diff(bins)/2; 
% y0(y0<nmin) = NaN;
% plot(x, movmean(100*y./y0,3, 'omitnan'), 'k+', 'Color', col, 'LineWidth',2); 
% quietaxes(gca, 'y');
leg = cleanLegend(gca, 'northeast', 'FontSize',fs-3);
leg.ItemTokenSize = leg.ItemTokenSize./2;

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
thisha = ha(4); axes(thisha); hold on; 
II = arrayfun(@(bi) find(isin((binvar), bins(bi:bi+1))), 1:length(bins)-1, 'Un',0);
% pp = find(cellfun(@(x) contains(x, 'PIC2'), proc));
% II = cellfun(@(ii) intersect(ii,pp), II,'Un',0);

x = bins(1:end-1) + diff(bins)/2;

var = swh; 
y = cellfun(@(ii) 100*sum(isnan(var(ii)))./length(ii),II);
y(cellfun(@length,II)<nmin) = NaN;
y = movmean(y,2, 'omitnan');
plot(x,y,'k:', 'MarkerSize',6, 'Color', [0.3 0.3 0.3]+0.2, 'LineWidth',2)

% 'suspect_beam_used suspect_less_than_nine_beams suspect_rain_likely suspect_pixel_used suspect_num_pt_avg suspect_karin_telem suspect_orbit_control suspect_sc_event_flag suspect_tvp_qual suspect_volumetric_corr degraded_beam_used degraded_large_attitude degraded_karin_ifft_overflow bad_karin_telem bad_very_large_attitude bad_outside_of_range degraded bad_not_usable'
flag_meanings = 'suspect_beam_used suspect_less_than_nine_beams suspect_rain_likely suspect_pixel_used suspect_num_pt_avg suspect_karin_telem suspect_orbit_control suspect_sc_event_flag suspect_tvp_qual suspect_volumetric_corr degraded_beam_used degraded_large_attitude degraded_karin_ifft_overflow bad_karin_telem bad_very_large_attitude bad_outside_of_range degraded bad_not_usable';
flag_meanings = strsplit(flag_meanings, ' ');
flag_masks    = [8          16          32         128         256         512        1024        2048        4096        8192      131072      262144      524288    16777216    33554432   536870912  1073741824  2147483648];
var = swh_qual; 
nn = find(var==0); var = log2(var); var(nn) = 0;
var = floor(var); 

cmap = subsetcmap(cmocean('solar'), 4);
cmap = cmap./1.3 + min((1-max(cmap)./1.3));
flgnms = {'bad', 'degraded', 'suspect'};

Y = NaN(length(flgnms), length(x));
for fi=1:length(flgnms)
    flgnm = flgnms{fi};
    ff = find(cellfun(@(x) contains(x, flgnm), flag_meanings)); ff = flag_masks(ff); ff = log2(ff);
    y = cellfun(@(ii) 100*sum(sum(var(ii)==ff')>0)./length(ii), II);
    Y(fi,:) = y; 
end
nn = find(~sum(isnan(Y),1));
Y = cumsum(Y); 
Y = movmean(Y,3,2,'omitnan');
for fi=length(flgnms):-1:1
    y = Y(fi,:);
    flgnm = flgnms{fi};
    patch([0.002 0.01 x(nn) max(x)+2],[0 y(nn(1)) y(nn) 0],'k',...
        'FaceAlpha', 0.9, 'FaceColor',  cmap(fi,:), 'EdgeColor', 'none',...
        'DisplayName', ['SWH ' flgnm])
end

var = rainrate; 
y = cellfun(@(ii) 100*sum(var(ii)>0)./length(ii),II); 
y(cellfun(@length,II)<nmin) = NaN; 
y = movmean(y,2, 'omitnan'); 
plot(x,y,'b.-', 'MarkerSize',20, 'DisplayName', 'rain', 'Color', [0.3 0.3 0.7])
leg = cleanLegend(gca, 'northeast', 'FontSize', fs-3);
leg.ItemTokenSize(1) =  leg.ItemTokenSize(1)/2;


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
thisha = ha(1); axes(thisha); hold on; 
bins = [-2.71:0.22:3.5]; binvar = disttocoast./1000; bins = 10.^bins; 
% yyaxis left; set(gca, 'YColor', 'k'); 
histogram(binvar(~isnan(hs)),bins, 'FaceColor', 'k', 'DisplayName', 'coastal buoy', 'EdgeColor', 'none')
y0 = histcounts(binvar(~isnan(hs)),bins);
histogram(binvar(~isnan(swh) & ~isnan(hs)), bins, 'FaceColor',  [0.3 0.75 0.9], 'DisplayName', 'SWOT KaRIn', 'EdgeColor', 'none')
y = histcounts(binvar(~isnan(swh) & ~isnan(hs)),bins);
% yyaxis right; col = [0.7 0.1 0.2]; set(gca, 'YColor', col); ylim([0 100])
% x = bins(1:end-1) + diff(bins)/2; 
% y0(y0<nmin) = NaN;
% plot(x, movmean(100*y./y0,3, 'omitnan'), 'k+', 'Color', col, 'LineWidth',2); 
% ylabel('\%', 'Interpreter', 'latex')
leg = cleanLegend(gca, 'northeast', 'FontSize',fs-3);
leg.ItemTokenSize = leg.ItemTokenSize./2;


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
thisha = ha(3); axes(thisha); hold on; 
bins = log10(bins);
II = arrayfun(@(bi) find(isin(log10(binvar), bins(bi:bi+1))), 1:length(bins)-1, 'Un',0);
% pp = find(cellfun(@(x) contains(x, 'PIC2'), proc));
% II = cellfun(@(ii) intersect(ii,pp), II,'Un',0);

x = bins(1:end-1) + diff(bins)/2; x = 10.^x; 

var = swh; 
y = cellfun(@(ii) 100*sum(isnan(var(ii)))./length(ii),II); y = movmean(y,2, 'omitnan'); 
y(cellfun(@length,II)<nmin) = NaN;
plot(x,y,'k:', 'MarkerSize',6, 'Color', [0.3 0.3 0.3]+0.2, 'LineWidth',2)

% 'suspect_beam_used suspect_less_than_nine_beams suspect_rain_likely suspect_pixel_used suspect_num_pt_avg suspect_karin_telem suspect_orbit_control suspect_sc_event_flag suspect_tvp_qual suspect_volumetric_corr degraded_beam_used degraded_large_attitude degraded_karin_ifft_overflow bad_karin_telem bad_very_large_attitude bad_outside_of_range degraded bad_not_usable'
flag_meanings = 'suspect_beam_used suspect_less_than_nine_beams suspect_rain_likely suspect_pixel_used suspect_num_pt_avg suspect_karin_telem suspect_orbit_control suspect_sc_event_flag suspect_tvp_qual suspect_volumetric_corr degraded_beam_used degraded_large_attitude degraded_karin_ifft_overflow bad_karin_telem bad_very_large_attitude bad_outside_of_range degraded bad_not_usable';
flag_meanings = strsplit(flag_meanings, ' ');
flag_masks    = [8          16          32         128         256         512        1024        2048        4096        8192      131072      262144      524288    16777216    33554432   536870912  1073741824  2147483648];
var = swh_qual; 
nn = find(var==0); var = log2(var); var(nn) = 0;
var = floor(var); 

cmap = subsetcmap(cmocean('solar'), 4);
cmap = cmap./1.3 + min((1-max(cmap)./1.3));
flgnms = {'bad', 'degraded', 'suspect'};

Y = NaN(length(flgnms), length(x));
for fi=1:length(flgnms)
    flgnm = flgnms{fi};
    ff = find(cellfun(@(x) contains(x, flgnm), flag_meanings)); ff = flag_masks(ff); ff = log2(ff);
    y = cellfun(@(ii) 100*sum(sum(var(ii)==ff')>0)./length(ii), II);
    Y(fi,:) = y; 
end
nn = find(~sum(isnan(Y),1));
Y = cumsum(Y); 
Y = movmean(Y,3,2,'omitnan');
for fi=length(flgnms):-1:1
    y = Y(fi,:);
    flgnm = flgnms{fi};
    patch([0.002 0.01 x(nn) max(x)+2],[0 y(nn(1)) y(nn) 0],'k',...
        'FaceAlpha', 0.9, 'FaceColor',  cmap(fi,:), 'EdgeColor', 'none',...
        'DisplayName', ['SWH ' flgnm])
end

var = rainrate; 
y = cellfun(@(ii) 100*sum(var(ii)>0)./length(ii),II); y = movmean(y,2, 'omitnan'); 
plot(x,y,'b.-', 'MarkerSize',20, 'DisplayName', 'rain', 'Color', [0.3 0.3 0.7])

leg = cleanLegend(gca, 'northeast', 'FontSize', fs-3);
leg.ItemTokenSize(1) =  leg.ItemTokenSize(1)/2;
AddLetters2Plots({ ha(1) ha(2) ha(3) ha(4)},...
     {'(a)', '(b)', '(c)', '(d)', '(e)', '(f)', '(g)'},...
    'BackgroundColor', 'none', 'Margin', 1,...
    'HShift', 0.01, 'VShift', 0.015, ...
    'FontName', 'Times', 'FontSize', 10, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top')
% savejpg(gcf, 'swot_vs_buoy_quality_vs_dist_and_hs', savepath, saveopt)

% savejpg(gcf, ['fig_datacoverage'], [base_path(1:12) 'desktop/'], 'on')

