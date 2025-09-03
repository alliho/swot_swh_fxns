base_path = '/Users/ajho/Documents/JPL/';
local_path = [base_path 'papers/swot_swh_calval/swot_swh_fxns/'];
addpath([local_path 'swh_fxns/matlab/'])
addpath([local_path 'fig_code/addfxns/'])
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

%% [save] JUST SCATTER | buoys and moorings | with nadir? no model 


mdl = correct_swotswh();
mdl = correct_swotswh('data', 'moorings');
% mdl = correct_swotswh('data', 'buoys');


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
scatopt = 0; clims = [10 70]
correctopt = 0; 
procopt = 1;
correctopt = 0; procopt = 1; showprocopt = 1; includepic2 = 1;
% correctopt = 1; procopt = 1; showprocopt = 1; includepic2 = 1;
addnadircenter = 1;
% -------------------------------------------------------------------------
dlabs = {'GNSS buoy', 'coastal buoy', 'GNSS buoy', 'coastal buoy'}
dds = [2 1 2 1];
fldnms = {'swh_nadir_altimeter','swh_nadir_altimeter', 'swh_karin_uncorr_unsmoothed', 'swh_karin_uncorr_unsmoothed'};
fldnms = {'swh_nadir_altimeter','swh_nadir_altimeter', 'swh_karin_uncorr_patch_unsmoothed', 'swh_karin_uncorr_patch_unsmoothed'};
fldnms = {'swh_nadir_altimeter','swh_nadir_altimeter', 'swh_karin_uncorr_patch', 'swh_karin_uncorr_patch'};
% fldnms = {'swh_nadir_altimeter','swh_nadir_altimeter', 'swh_karin_uncorr', 'swh_karin_uncorr'};
fldlbs = {'SWOT nadir','SWOT nadir', 'SWOT KaRIn', 'SWOT KaRIn'};



qualflag = 1;
uncflag = 1;
% -------------------------------------------------------------------------
figure(679174); clf; 
setfigsize(gcf, [662   627])
% tight_subplot(Nh, Nw, [gap_h gap_w], [lower upper], [left right]) 
% oldha = tight_subplot(2, 4, [0.025 0.05], [0.11 0.05], [0.05 0.05]);
% [ha] = reorg_tightsubplot(oldha, {[1 2],[4 5]+1, [3 7], [4 8
% ha = tight_subplot(3, 2, [0.025 0.05], [0.12 0.025], [0.033 0.033]);
ha = tight_subplot(2, 2, [0.02 0.02], [0.08 0.015], [0.07 0.07]);
% [ha] = reorg_tightsubplot(oldha, {[1 2],[4 5], [3 6]});

% -------------------------------------------------------------------------

setaxes(ha(:),'XLim', [0 6.2]);
setaxes(ha(:),'YLim', [0 6.2]);
setaxes(ha(1:end),'XLim', [0 7.5]);
setaxes(ha(1:end),'YLim', [0 7.5]);
setaxes(ha(:),'CLim', [10 60]);
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
    q = cell2mat(cellfun(@(x) x(:)', {data(dd).swh_karin_qual}, 'Un', 0));
    x = cell2mat(cellfun(@(x) x(:)', {data(dd).hs}, 'Un', 0));
    y = cell2mat(cellfun(@(x) x(:)', {data(dd).(fldnm)}, 'Un', 0)); y(y<0) = NaN;
    z = cell2mat(cellfun(@(x) x(:)', {data(dd).cross_track_distance}, 'Un', 0)); z = abs(z);

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
    

    if addnadircenter & contains(fldnm, 'nadir') & dd~=2
        nn = find(z<=2 & ~isnan(x) & ~isnan(y)); 
        % nn = intersect(nn , find(cellfun(@(x) strcmp(x, 'PIC2'), proc))); 
        x_nadir = x(nn); 
        y_nadir = y(nn); 
    else
        x_nadir = NaN; y_nadir = NaN; 
    end


    nn = cellfun(@(fldnm) isnan(cell2mat(cellfun(@(x) x(:), {data(dd).(fldnm)}, 'Un', 0))), fldnms, 'Un', 0);
    nn = cell2mat(nn); nn = sum(nn,2); nn = find(~nn); % only keep data that overlaps in coverage
    A = {x,y,z,hs, swh_nadir, swh_model, xtrk, proc, depth, disttocoast,swh_msk,swh_raw,swh_qual, swh_proc};
    A = cellfun(@(x) x(nn), A, 'Un', 0);
    [x,y,z,hs, swh_nadir, swh_model, xtrk, proc, depth, disttocoast,swh_msk,swh_raw,swh_qual, swh_proc] = A{:};

    

    % only look inside tracks
    nn = find(isin(xtrk,[10 60]));
    A = {x,y,z,hs, swh_nadir, swh_model, xtrk, proc, depth, disttocoast,swh_msk,swh_raw,swh_qual, swh_proc};
    A = cellfun(@(x) x(nn), A, 'Un', 0);
    [x,y,z,hs, swh_nadir, swh_model, xtrk, proc, depth, disttocoast,swh_msk,swh_raw,swh_qual, swh_proc] = A{:};
    
    
    % remove bad data (zero values)
    nn = find(swh_raw~=0); 
    A = {x,y,z,hs, swh_nadir, swh_model, xtrk, proc, depth, disttocoast,swh_msk,swh_raw,swh_qual, swh_proc};
    A = cellfun(@(x) x(nn), A, 'Un', 0);
    [x,y,z,hs, swh_nadir, swh_model, xtrk, proc, depth, disttocoast,swh_msk,swh_raw,swh_qual, swh_proc] = A{:};
    
    % remove outliars 
    nn = find(abs(swh_proc-hs)./hs < 1); 
    A = {x,y,z,hs, swh_nadir, swh_model, xtrk, proc, depth, disttocoast,swh_msk,swh_raw,swh_qual, swh_proc};
    A = cellfun(@(x) x(nn), A, 'Un', 0);
    [x,y,z,hs, swh_nadir, swh_model, xtrk, proc, depth, disttocoast,swh_msk,swh_raw,swh_qual, swh_proc] = A{:};
    
    % -------------------------------------------------------------------------
    
    if contains(fldnm,'swh_karin') & correctopt
        corr = mdl(y, z); corr(isinf(corr)) = 0;
        nn = find(cellfun(@(x) strcmp(x, 'PIC2'), proc)); 
        corr(nn) = 0; 
        corr(isinf(corr)) = NaN;
        y = y - corr;
    else
        corr = zeros(size(x));
    end

    

    % A = {xdata,ydata,zdata,pdata}; 
    % nn = getnnn(xdata,ydata,zdata);
    % A = cellfun(@(x) x(nn), A, 'Un', 0);
    % [xdata,ydata,zdata,pdata] = A{:};

    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    setaxes(ha(hi),'YLabel', [ylb ' SWH [m]']);
    setaxes(ha(hi),'XLabel', [xlb ' $H_s$ [m]'] );
    clear sts
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    thisha = ha(hi); axes(thisha); hold on;

    
    if scatopt
        scatter(x, y,15, z, 'filled')
        if hi==length(ha)              
           c = fixedcolorbar(gca, 'Location', 'eastoutside');
        end
    elseif procopt & ~scatopt & contains(fldnm, 'karin')
        nn = find(cellfun(@(x) strcmp(x, 'PIC0') | strcmp(x, 'PGC0'), proc)); 
        procname = 'P..0';
        nn = ~cellfun(@(x) ~isempty(regexp(x,procname)), proc); nn = find(nn==0); 
        labstr = 'Version C'
        scatter(x(nn), y(nn),6, 'k', 'filled', 'DisplayName', labstr, 'MarkerFaceAlpha', 0.6)

        procname = 'PIC2';
        nn = ~cellfun(@(x) ~isempty(regexp(x,procname)), proc); nn = find(nn==0); 

        if ~isempty(nn)
            if showprocopt 
                scatter(x(nn), y(nn) ,7, [0.6 0 0],'filled', '^', 'DisplayName', 'Version D')
            elseif includepic2
                scatter(x(nn), y(nn),6, 'k', 'filled', 'MarkerFaceAlpha', 0.6)
            end
        end
        
    else
        scatter(x, y, 6, 'k', 'filled', 'MarkerFaceAlpha', 0.5)
    end


    %%% ADD STATS LABELS
    dx = - 0.12; dy = 0.015;
    if ~procopt | ~contains(fldnm, 'karin')
        xdata = x; ydata = y; 
        sts.cc = ccfxn(xdata,ydata);
        sts.rmse = rmsefxn(xdata,ydata);
        sts.crmse = crmsefxn(xdata,ydata);
        sts.bias = biasfxn(xdata,ydata);
        sts.SI = sifxn(xdata,ydata);
        x0 = sum(thisha.Position([1 3])) + dx; y0 = thisha.Position(2) + dy; 
        textbypos(x0,y0,...
            {   ['n$=$' num2str(sum(~isnan(xdata) & ~isnan(ydata)))], ...
                ['bias$=$' num2str(round(nanmean([sts.bias]),2)) 'm'], ...
                ['cRMSE$=$' num2str(round(nanmean([sts.crmse   ]),2)) 'm'],...
                ['CC$=$' num2str(round(nanmean([sts.cc]),2))]}, ...
            'interpreter', 'latex', 'VerticalAlignment','bottom');

    elseif procopt & contains(fldnm, 'karin')
        if ~showprocopt & correctopt & includepic2
            % nn = find(cellfun(@(x) strcmp(x, 'PIC0') | strcmp(x, 'PGC0') |  strcmp(x, 'PIC2'), proc)); 
            procname = 'P...';
            nn = ~cellfun(@(x) ~isempty(regexp(x,procname)), proc); nn = find(nn==0); 
        else
            procname = 'P..0';
            nn = ~cellfun(@(x) ~isempty(regexp(x,procname)), proc); nn = find(nn==0); 
        end
        xdata = x(nn); ydata = y(nn);

        sts.cc = ccfxn(xdata,ydata);
        sts.rmse = rmsefxn(xdata,ydata);
        sts.crmse = crmsefxn(xdata,ydata);
        sts.bias = biasfxn(xdata,ydata);
        sts.SI = sifxn(xdata,ydata);
        x0 = sum(thisha.Position([1 3])) + dx; y0 = thisha.Position(2) + dy; 
        textbypos(x0,y0,...
            {   ['n$=$' num2str(sum(~isnan(xdata) & ~isnan(ydata)))], ...
                ['bias$=$' num2str(round(nanmean([sts.bias]),2)) 'm'], ...
                ['cRMSE$=$' num2str(round(nanmean([sts.crmse   ]),2)) 'm'],...
                ['CC$=$' num2str(round(nanmean([sts.cc]),2))]}, ...
            'interpreter', 'latex', 'VerticalAlignment','bottom');

        nn = find(cellfun(@(x) strcmp(x, 'PIC2'), proc)); 
        if ~isempty(nn) & showprocopt & contains(fldnm, 'karin')
            xdata = x(nn); ydata = y(nn);
            sts.cc = ccfxn(xdata,ydata);
            sts.rmse = rmsefxn(xdata,ydata);
            sts.crmse = crmsefxn(xdata,ydata);
            sts.bias = biasfxn(xdata,ydata);
            sts.SI = sifxn(xdata,ydata);
            x0 = sum(thisha.Position([1 3])) + dx - 0.12; y0 = thisha.Position(2) + dy; 
            if showprocopt & contains(fldnm, 'karin')
            textbypos(x0,y0,...
                {   ['n$=$' num2str(sum(~isnan(xdata) & ~isnan(ydata)))], ...
                    ['bias$=$' num2str(round(nanmean([sts.bias]),2)) 'm'], ...
                    ['cRMSE$=$' num2str(round(nanmean([sts.crmse   ]),2)) 'm'],...
                    ['CC$=$' num2str(round(nanmean([sts.cc]),2))]}, ...
                'interpreter', 'latex', 'VerticalAlignment','bottom', 'Color', [0.6 0 0]);
            end
        end




        
    end

    if ~isnan(x_nadir)
        col = [0.11 0.28 0.8]+[0.3 0.3 0.2];
        scatter(x_nadir, y_nadir ,11, col,'filled', '^', 'DisplayName', '\Deltax < 2km')
        xdata = x_nadir; ydata = y_nadir;
        sts.cc = ccfxn(xdata,ydata);
        sts.rmse = rmsefxn(xdata,ydata);
        sts.crmse = crmsefxn(xdata,ydata);
        sts.bias = biasfxn(xdata,ydata);
        sts.SI = sifxn(xdata,ydata);
        x0 = sum(thisha.Position([1 3])) + dx - 0.12; y0 = thisha.Position(2) + dy; 
        
        textbypos(x0,y0,...
            {   ['n$=$' num2str(sum(~isnan(xdata) & ~isnan(ydata)))], ...
                ['bias$=$' num2str(round(nanmean([sts.bias]),2)) 'm'], ...
                ['cRMSE$=$' num2str(round(nanmean([sts.crmse   ]),2)) 'm'],...
                ['CC$=$' num2str(round(nanmean([sts.cc]),2))]}, ...
            'interpreter', 'latex', 'VerticalAlignment','bottom', 'Color', col);


    end

    if (showprocopt & contains(fldnm, 'karin')) | ~isnan(x_nadir)
        cleanLegend(gca, 'northeast', 'FontSize',9, 'NumColumns',2)
    end
    
end


quietaxes(ha(setdiff([1:4],[3 4])), 'x')
quietaxes(ha(setdiff([1:4],[1:2:4])), 'y')

AddLetters2Plots({ ha(1) ha(2) ha(3) ha(4)},...
     {'(a)', '(b)', '(c)', '(d)', '(e)', '(f)', '(g)'},...
    'BackgroundColor', 'w', 'Margin', 1,...
    'HShift', 0.01, 'VShift', 0.015, ...
    'FontName', 'Times', 'FontSize', 10, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top')

optstr = '';
if scatopt; optstr = [optstr '_scatter']; end
if correctopt; optstr = [optstr '_corrected']; end
if procopt; optstr = [optstr '_wproc']; end
if addnadircenter; optstr = [optstr '_wcenternadir']; end


