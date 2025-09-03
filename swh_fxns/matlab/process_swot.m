function [var] = process_swot(swot, fldnm, varargin)
% process_swot clean up swot
% [fld] = process_swot(swot, fldnm, varargin)
% 
%    modopt = 0; smoothopt = 1; ri = 2; corropt = 1; patchopt = 1; mskopt = 1; 
%    varargin = {'patch', patchopt, 'model', mdl, 'smooth', smoothopt, 'ri', ri, 'correct', corropt, 'mask', mskopt};
%    [Z] = process_swot(swot, fldnm, varargin{:});

%% set inherent variables
dr = 2; % grid spacing 2 km  

%% set default variables
minflag = 1; mindepth = -1; ri = 3; Lavg = ri*dr;
swhopt = 0; sshopt = 0; wspopt = 0; 
modopt = 0; smoothopt = 0; rainopt = 0;
patchopt = 0; corropt = 0; mskopt = 0; 
dspkopt = 0;
ptchsz = 3;

if contains(fldnm, 'ssh'); sshopt = 1; end
if contains(fldnm, 'swh'); swhopt = 1; end
if contains(fldnm, 'wind_speed'); wspopt = 1; end

%% parse variable inputs
p = inputParser; p.KeepUnmatched = true;
addParameter(p, 'swh', swhopt);
addParameter(p, 'qcmodel', modopt);
addParameter(p, 'model', NaN);
addParameter(p, 'rain', rainopt);
addParameter(p, 'wind_speed', wspopt);
addParameter(p, 'ssha', sshopt);
addParameter(p, 'ssh', sshopt);
addParameter(p, 'smooth', smoothopt);
addParameter(p, 'mindepth', mindepth);
addParameter(p, 'ri', ri);
addParameter(p, 'Lavg', Lavg);
addParameter(p, 'patch', patchopt);
addParameter(p, 'correct', corropt);
addParameter(p, 'mask', mskopt);
addParameter(p, 'despike', dspkopt);
addParameter(p, 'patchsize', ptchsz);
addParameter(p, 'quality', 'good');
parse(p, varargin{:});
if p.Results.ssha~=sshopt
    sshopt = p.Results.ssha;
elseif p.Results.ssh~=sshopt
    sshopt = p.Results.ssh;
end
wspopt = p.Results.wind_speed;
rainopt = p.Results.rain;
swhopt = p.Results.swh;
mdl = p.Results.model;
patchopt = p.Results.patch;
modopt = p.Results.qcmodel;
corropt = p.Results.correct;
mskopt = p.Results.mask;
ri = p.Results.ri;
Lavg = p.Results.Lavg;
mindepth = p.Results.mindepth;
smoothopt = p.Results.smooth;
dspkopt = p.Results.despike;
ptchsz = p.Results.patchsize;
qualflag = p.Results.quality;

if isa(mdl, 'NaN'); corropt = 0; end
if patchopt==2; dspkopt=1; end
if mskopt>0; rainopt=1; end

if ~any(strcmp(varargin, 'Lavg')) & any(strcmp(varargin, 'ri'))
    Lavg = ri*dr;
elseif any(strcmp(varargin, 'Lavg')) & ~any(strcmp(varargin, 'ri'))
    ri = Lavg./dr; ri = ceil(ri);
end



%% get swot label
sfldnm = strrep(fldnm, '_qc_avg', '');
sfldnm = strrep(sfldnm, '_qc', '');
sfldnm = strrep(sfldnm, '_cor', '');
sfldnm = strrep(sfldnm, '_avg', '');
sfldnm = strrep(sfldnm, '_proc', '');

if contains(fldnm, 'height_cor_xover'); sfldnm = fldnm; end

%% determine field and correct accordingly



if strcmp(sfldnm, 'time')
    var = swot.(sfldnm); % ASSIGN
    return
end


% if contains(fldnm, 'ssha') | contains(fldnm, 'ssh')
%     var = swot.(sfldnm) + swot.height_cor_xover; 
% else
%     var = swot.(sfldnm);
% end

var = swot.(sfldnm); 
var(isinf(var)) = NaN;
rawvar = var;


if mskopt
    msk = ones(size(var));
    msk(isnan(var)) = NaN; 
    
    if isfield(swot, [sfldnm '_qual'])
        quvar = swot.([sfldnm '_qual']);
        quinfo = ncinfo([swot.fpath swot.fname], [sfldnm '_qual']);
        flag_meanings = quinfo.Attributes(4).Value; flag_meanings = strsplit(flag_meanings);
        flag_masks = quinfo.Attributes(5).Value;
        
        if strcmp(qualflag, 'good')
            msk(quvar > 0) = NaN;
        else
            % flag_heirarchy = {'suspect', 'degraded', 'bad'};
            qq = find(contains(flag_meanings, qualflag));
            maxflag = max(flag_masks(qq));
            msk(quvar > maxflag) = NaN;
        end
        
        % if contains(sfldnm, 'swh_karin')
        %     msk(isin(log2(quvar), [11 11.1], 'inclusive')) = 1; % gets rid of the calibration patches
        %     msk(isnan(swot.([sfldnm '_uncert']))) = NaN;
        % elseif contains(sfldnm, 'ssh') & contains(sfldnm, 'karin')
        %     msk(isin(log2(quvar), [11 11.1], 'inclusive')) = 1; % gets rid of the calibration patches
        %     % msk(isnan(swot.([sfldnm '_uncert']))) = NaN;
        % end
    
        if mskopt>=2
            msk(isin(log2(quvar), [11 11.2], 'inclusive')) = 1; % gets rid of the calibration patches
            msk(isnan(swot.([sfldnm '_uncert']))) = NaN;
            if contains(sfldnm, 'swh_karin') & mskopt>=3  
                msk(isin(log2(quvar), [7 7.2], 'inclusive')) = 1;
                msk(isin(log2(quvar), [5 7], 'inclusive')) = 1; % gets rid of the rain patches..?
            end
            if contains(sfldnm, 'wind_speed_karin') & mskopt>=3 
                % msk(isin(log2(quvar), [7 7.2], 'inclusive')) = 1;
                msk(isin(log2(quvar), [5 7], 'inclusive')) = 1; % gets rid of the rain patches..?
            end
        end
    end
    

%     if contains(sfldnm, 'wind_speed')
%         [msk] = qc_swot(swot,'wind_speed',1, 'swh',0, 'model', 0,...
%                         'smooth', smoothopt, 'ri', ri, 'mindepth', mindepth);
%     elseif contains(sfldnm, 'ssh')   
%         [msk] = qc_swot(swot,'ssh',1,        'swh',0, 'model', 0,...
%                         'smooth', smoothopt, 'ri', ri, 'mindepth', mindepth);
%     else                              
%         [msk] = qc_swot(swot, 'model', modopt,...
%                         'smooth', smoothopt, 'ri', ri, 'mindepth', mindepth);
%     end 
%     if mskopt==2
%         msk(isin(log2(swot.swh_karin_qual), [11 11.2], 'inclusive')) = 1; % gets rid of the calibration patches
%         if contains(sfldnm, 'swh'); 
%             msk(isin(log2(swot.swh_karin_qual), [7 7.2], 'inclusive')) = 1;
%         end
%     end


    
    msk(swot.depth_or_elevation > mindepth) = NaN;
    if rainopt
        rainmsk = ones(size(var));
        % rainmsk(swot.rain_rate > 1) = 0;
        rainmsk(swot.rain_rate > 0) = 0;
        rainmsk(swot.rain_flag > 0) = 0;
        % rainmsk = filter2(ones(3)./9, rainmsk); rainmsk(rainmsk~=1) = NaN;
        % % rainmsk = filter2(ones(2)./4, rainmsk); 
        % rainmsk = round(rainmsk, 1); rainmsk(rainmsk<0.33) = NaN; rainmsk(~isnan(rainmsk)) = 1;
        rainmsk(rainmsk==0) = NaN;
        msk = msk.*rainmsk;
    end

    if contains(fldnm, 'swh_karin'); 
        msk(var<=0) = NaN;
    end
    var = var.*msk;
end

if dspkopt
    pi = 2; conn = 4; 
    msk = ones(size(var));
    msk(isnan(var)) = NaN; 
    msk_anomalies = makeanomalymask(var, pi, conn);
    var = var.*msk_anomalies;
    var = fillmissing2(var, 'linear');
    var = var.*msk;
    % if contains(fldnm, 'swh_karin'); var(var<0) = NaN; end
end


if patchopt==1
    var = fillmissing2(var,"movmedian", ri);
elseif patchopt==2
    % anomaly patching!!!
    pi = ptchsz; conn = 4; 
    if contains(fldnm, 'swh'); 
        msk_anomalies = makeanomalymask(var, pi, conn);
        var = var.*msk_anomalies;
    end
    var = patchswot(var, pi-1, conn);
    var = patchswot(var, pi+1, conn);
end



if smoothopt    
    msk = ones(size(var)); msk(isnan(var)) = NaN;

    
    
    % msk = ~isnan(var); 
    % pi = (2*ri)^2; conn = 4; 
    % msk = abs(bwareaopen(abs(msk-1),pi,conn)-1); % create mask that allows NaNs that are small
    % msk = abs(bwareaopen(abs(msk),pi+1,conn)); % and remove data that is small
    % msk = double(msk); msk(msk~=1) = NaN; 
    % msk_anomalies = ones(size(var));

    % var = fillmissing2(var.*msk_anomalies, 'linear'); % linear fill but skip over anomalous single points? 
    % H = fspecial('gaussian',3,ri/7);
    % var = filter2(H,var);
    % H = fspecial('gaussian',3,ri/3.5);
    % var = filter2(H,var);
    % var = var.*msk; 


    % %%% OLD? 
    % % H = ones(ri,ri)/ri^2; 
    % % var = conv2(var,H,'same'); 
    % H = fspecial('gaussian',ri,ri/1.5);
    % % H = ones(ri,ri)/ri^2; 
    % % H = hann(ri)*hann(ri)';
    % var = filter2(H,var);
    % var = var.*msk;


    %%% NEW WHILE PRESCRIBING FWHM BASED ON LAVG 2025/07/20
    N = floor((Lavg/2 + 1)./2)*3; 
    N = floor(N/2)*2+1;
    FWHM = Lavg;
    sigma = FWHM / (2 * sqrt(2 * log(2)));
    ri = sigma/2; 
    H = fspecial('gaussian',N,ri);
    var = filter2(H,var);
    var = var.*msk;


end

if patchopt==1
    % var = fillmissing2(var,"movmedian", ri);
    var = fillmissing2(var,"movmean", ri);
end

if smoothopt
    % remove center that would have been filled in by smoothing and patch fill.
    crosstracklims = [10 60] + [-1 1]*(ri*2000)/2;
    msk = NaN(size(var)); 
    msk(isin(abs(swot.cross_track_distance./1000), crosstracklims, 'inclusive') & swot.depth_or_elevation<mindepth) = 1; 
    var = var.*msk; 
end


if corropt & contains(fldnm, 'swh_karin') & ~contains(swot.processing, 'PIC2')
    % mdl = correct_swotswh();
    msk = ones(size(var)); 
    msk(isnan(var) | var < 0.05) = NaN;
    corr = mdl(var, abs(swot.cross_track_distance./1000)); corr(isinf(corr)) = 0; 
    var = var - corr; 
    var = var.*msk;
    var(var<0) = NaN;

    if patchopt==2
        % anomaly patching!!!
        pi = 3; conn = 4; 
        if contains(fldnm, 'swh'); 
            msk_anomalies = makeanomalymask(var, pi, conn);
            var = var.*msk_anomalies;
        end
        var = patchswot(var, pi, conn);
    end



end



end

function var = patchswot(var, pi, conn)
    cc = find(sum(isnan(var),2)==size(var,2));

    msk = ~isnan(var); 
    msk(cc,:) = 1;
    msk = abs(bwareaopen(abs(msk-1),pi,conn)-1); % create mask that allows NaNs that are small
    msk = abs(bwareaopen(abs(msk),pi+1,conn)); % and remove data that is small
    msk = double(msk); msk(msk~=1) = NaN;     
    var = fillmissing2(var, 'linear'); % linear fill but skip over anomalous single points? 
    msk(cc,:) = NaN;
    var = var.*msk; 

end

function msk_anomalies = makeanomalymask(var, pi, conn)
    msk_anomalies = ones(size(var));
    var_large = filt2(var,2000,15*1000, 'lp');
    var_small = filt2(var,2000,15*1000, 'hp');
    msk_anomalies = double(msk_anomalies); 
    msk_anomalies((var_small./var_large) > 0.2) = NaN;

    
end

