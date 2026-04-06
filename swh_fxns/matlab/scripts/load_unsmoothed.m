fpath = '/Users/ajho/Downloads/';
fname =     'SWOT_L2_LR_SSH_Unsmoothed_009_433_20240119T135934_20240119T145101_PGD0_02.nc';
swot = load_any_nc([fpath fname]); 


lonlims = 167.481015 + [-1 1].*10; lonlims = lonlims - 360;
latlims = 9.431479 + [-1 1].*10;

figure(802); clf; hold on; 
grps = {'left_', 'right_'};
for grp = grps
    grp = grp{1};
    X = swot.([grp 'longitude']); Y = swot.([grp 'latitude']);
    Z = swot.([grp 'ssha_karin_2']);

    % add cross-track correciton
    Z = Z + swot.([grp 'height_cor_xover']);
    
    % mask out bad data
    msk = ones(size(Z)); 
    msk(swot.([grp 'ssha_karin_2_qual']) > 0) = NaN;

    mskvar = deg2km(mag(swot.([grp 'latitude_uncert']), swot.([grp 'longitude_uncert']))).*1000;    
    msk(mskvar > 2.*repmat(nanmedian(mskvar,2), [1, size(mskvar,2)])) = NaN; 
    % i use the lat/lon uncertainty to mask out bad pixels that for some reason aren't flagged in the _qual fields
    % here i'm masking out any pixels that have more than 2x the along-track median uncertainty
    % slightly arbitrary so can adjust as needed

    Z = Z.*msk;
    ii = find(isin(Y(30,:), latlims)); 
    X = X(:,ii); Y = Y(:,ii); Z = Z(:,ii); 
    pcolor(X,Y,Z); shading flat;
end
caxis([-1 1].*0.4); 
colormap(redblue)