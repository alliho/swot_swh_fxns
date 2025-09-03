base_path = '/Users/ajho/Documents/JPL/';
local_path = [base_path 'papers/swot_swh_calval/swot_swh_fxns/'];
addpath([local_path 'swh_fxns/matlab/'])
addpath([local_path 'fig_code/addfxns/'])
swot_fpath = [base_path '/data/SWOT/onedayrepeat/'];
load('/Users/ajho/Documents/myrepos/supportingdata/nicejet.mat');
mdl = correct_swotswh();

%% LOAD SWOT
% SWOT Expert Data has been downloaded to local folder |swot_fpath| via earthdata.nasa.gov
latlims = [31 38]; lonlims = [-129 -123]; daterange = [datenum(2023,4,18) datenum(2023,5,2)];
SWOT = load_swot(swot_fpath,lonlims, latlims, daterange);

%% PROCESS SWH DATA
si = 1; 
swot = SWOT(si); 

% -------------------------------------------------------------------------
figure(1); clf; 
setfigsize(gcf, [1495         335])
% tight_subplot(Nh, Nw, [gap_h gap_w], [lower upper], [left right]) 
ha = tight_subplot(1, 6, [0.01 0.008], [0.09 0.1], [0.03 0.01]);
% [ha] = reorg_tightsubplot(oldha, {[1],[2:3]});
% -------------------------------------------------------------------------
format_fig4print(ha);
setaxes(ha, 'DataAspectRatio', [1 1 1]); 
setaxes(ha, 'XLim', lonlims + [-1 1].*(-0.1)); 
setaxes(ha, 'YLim', latlims + [-1 1].*(-0.1)); 
setaxes(ha, 'CLim', [1.6 4.4]); 
setaxes(ha, 'CLim', [1.6 3.4]); 
% setaxes(ha, 'CLim', [2 4.6]); 
% setaxes(ha, 'CLim', [1 2]); 
% setaxes(ha, 'CLim', [2.5 4]); 
setaxes(ha, 'Colormap', nicejet); 
linkaxes(ha)
format_mapticks(ha, 'xy');
quietaxes(ha(2:end), 'y'); 

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
sgtitle({datestr(swot.t0), ''}, 'FontName','times', 'FontWeight', 'bold');
fs = 11; pi =  5;
quflag = 'good'; mskopt = 1;
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
X = swot.longitude - 360; Y = swot.latitude; fldnm = 'swh_karin';
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
thisha = ha(1); axes(thisha); hold on; 
Z = swot.(fldnm);
pcolor(X,Y,Z); shading flat;
optstr = {'swh\_karin'};
text(lonlims(1)+0.2, latlims(2) - 0.2, optstr, 'FontName','courier', 'FontSize',fs, 'VerticalAlignment','top')
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
thisha = ha(2); axes(thisha); hold on; 
Z = process_swot(swot, fldnm, 'mask',mskopt, 'quality', quflag);
pcolor(X,Y,Z); shading flat;
optstr = [optstr, 'GOOD QA'];
text(lonlims(1)+0.2, latlims(2) - 0.2, optstr, 'FontName','courier', 'FontSize',fs, 'VerticalAlignment','top')
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
thisha = ha(3); axes(thisha); hold on; 
Z = process_swot(swot, fldnm, 'mask',mskopt, 'quality', quflag, 'despike', 1);
pcolor(X,Y,Z); shading flat;
optstr = [optstr, 'de-spiked'];
text(lonlims(1)+0.2, latlims(2) - 0.2,optstr, 'FontName','courier', 'FontSize',fs, 'VerticalAlignment','top')
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
thisha = ha(4); axes(thisha); hold on; 
Z = process_swot(swot, fldnm, 'mask',1, 'quality', quflag, 'despike', 1, 'patch',2, 'patchsize',pi);
pcolor(X,Y,Z); shading flat;
optstr = [optstr, 'patched'];
text(lonlims(1)+0.2, latlims(2) - 0.2, optstr, 'FontName','courier', 'FontSize',fs, 'VerticalAlignment','top')

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
thisha = ha(5); axes(thisha); hold on; 
Z = process_swot(swot, fldnm, 'mask',mskopt, 'quality', quflag, 'despike', 1, 'patch',2, 'patchsize',pi, 'correct', 1, 'model', mdl);
pcolor(X,Y,Z); shading flat;
optstr = [optstr, 'empirically corrected'];
text(lonlims(1)+0.2, latlims(2) - 0.2, optstr, 'FontName','courier', 'FontSize',fs, 'VerticalAlignment','top')
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
thisha = ha(6); axes(thisha); hold on; 
Z = process_swot(swot, fldnm, 'mask',mskopt, 'quality', quflag, 'despike', 1, 'patch',2, 'patchsize',pi, 'correct', 1, 'model', mdl, 'smooth', 1, 'Lavg',5);
pcolor(X,Y,Z); shading flat;
optstr = [optstr, '5 km Gaussian filter'];
text(lonlims(1)+0.2, latlims(2) - 0.2, optstr, 'FontName','courier', 'FontSize',fs, 'VerticalAlignment','top')




savejpg(gcf, ['example_processing_' datestr(swot.t0, 'yyyymmddTHHMMSS')], '/Users/ajho/Documents/JPL/figs/', 'time')