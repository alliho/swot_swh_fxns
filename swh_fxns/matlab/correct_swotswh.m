function [mdl] = correct_swotswh(varargin)

% model syntax = mdl(y,z) = mdl(swh_karin, cross_track_distance)

%%
p = inputParser;
addParameter(p, 'data', []);
addParameter(p, 'processing', {'PGC0', 'PIC0'});
parse(p, varargin{:});
datatype = p.Results.data;
processingalgs = p.Results.processing;

%% load data
fpath = '/Users/ajho/Documents/JPL/data/processed/';
% tmp = load([fpath 'buoys_SWOT.mat']); data(1) = tmp.data;
% tmp = load([fpath 'MO_SWOT.mat']); data(2) = tmp.data;

tmp = load([fpath 'buoys_SWOT_w-wind.mat']); data(1) = tmp.data;
tmp = load([fpath 'MO_SWOT_w-wind.mat']); data(2) = tmp.data;
% tmp = load([fpath 'buoys_SWOT_onedayrepeat_w-wind.mat']); data(3) = tmp.data;



%% choose data
if strcmp(datatype, 'moorings')
    dd = 2;
elseif strcmp(datatype, 'buoys')   
    dd = 1;
else
    dd = 1:2;
end

%% correction
% rawfldnm = 'swh_karin_uncorr_unsmoothed';
rawfldnm = 'swh_karin_uncorr_patch';
% rawfldnm = 'swh_karin_uncorr_patch_unsmoothed';
lims = [1 99];

x = cell2mat(cellfun(@(x) x(:)', {data(dd).hs}, 'Un', 0));
y = cell2mat(cellfun(@(x) x(:)', {data(dd).(rawfldnm)}, 'Un', 0));
z = cell2mat(cellfun(@(x) x(:)', {data(dd).cross_track_distance}, 'Un', 0)); z = abs(z);
p = cellfun(@(x) x(:)', {data(dd).processing}, 'Un', 0); p = [p{:}]';
nn = find(~isnan(x) & ~isnan(y) & ~isnan(z)); x = x(nn); y = y(nn); z = z(nn); p = p(nn);
nn = find(isin(y-x, prctile(y-x, lims))); x = x(nn); y = y(nn); z = z(nn); p = p(nn);

hs = cell2mat(cellfun(@(x) x(:)', {data(dd).hs}, 'Un', 0));
swh = cell2mat(cellfun(@(x) x(:)', {data(dd).(rawfldnm)}, 'Un', 0));
xtrk = cell2mat(cellfun(@(x) abs(x(:))', {data(dd).cross_track_distance}, 'Un', 0)); 
proc = cellfun(@(x) x(:)', {data(dd).processing}, 'Un', 0); proc = [proc{:}]';

nn = find(isin(abs(swh-hs), prctile(abs(swh-hs), lims)) & ~isnan(hs) & ~isnan(swh)); 
A = {hs, swh, xtrk, proc};
A = cellfun(@(x) x(nn), A, 'Un', 0);
[hs, swh, xtrk,proc] = A{:};


%% choose processing algorithm
% processingalgs

nn = find(cellfun(@(x) strcmp(x, 'PGC0') | strcmp(x, 'PIC0'), proc ));
A = {hs, swh, xtrk, proc};
A = cellfun(@(x) x(nn), A, 'Un', 0);
[hs, swh, xtrk,proc] = A{:};

opts =  optimset('display','off');

%% uniformly subsample in SWH and XTRCK
x = swh'; y = xtrk';
prckeep = 0.75; 
xbins = [0:0.25:10];
ybins = [0:5:90];

Nbx = numel(xbins) - 1;
Nby = numel(ybins) - 1;

xb = discretize(x, xbins);
yb = discretize(y, ybins);

Nb = accumarray([xb yb], 1, [Nbx Nby]);

% Combine 2D bin indices into a single group ID
binID = (xb - 1) * (numel(ybins) - 1) + yb;


% Group IDs — keep NaNs; they'll be ignored in grouping
[~, ~, binGroups] = unique(binID, 'rows', 'stable');

% Preallocate
Nkeep = max([prctile(Nb(:),95).*prckeep,2]);
Nkeep = max([prctile(Nb(:),95).*prckeep]); Nkeep = round(Nkeep);
idx_sub = [];

% Loop over unique bin groups (ignoring NaN bin IDs)
uniqueBins = unique(binID(~isnan(binID)));
for i = 1:numel(uniqueBins)
    bb = find(binID == uniqueBins(i));
    n = min(numel(bb), Nkeep);
    idx_sub = [idx_sub; randsample(bb, n)];
end

% Final sampled subset
% x_sub = x(idx_sub);
% y_sub = y(idx_sub);

nn = idx_sub;
A = {hs, swh, xtrk, proc};
A = cellfun(@(x) x(nn), A, 'Un', 0);
[hs, swh, xtrk, proc] = A{:};



%% model
% mdl = fit([y; z]', y'-x','poly22');
% mdl = fit([y; z]', y'-x','poly22');

% ft = fittype('(b1*x + b2)*((y-c1)^2 + c2) + a',...
%     'dependent',{'z'},'independent',{'x', 'y'},...
%     'coefficients',{'a', 'b1', 'b2', 'c1', 'c2'});
% options = fitoptions('Method', 'NonlinearLeastSquares', ...
%                          'Lower',[-Inf -Inf 0 10 -Inf], 'Upper',[Inf 0 Inf 60 Inf]);
% mdl = fit([y; z]', y'-x',ft)

% 
% ft = fittype('p00 + p10*x + p01*y + p20*x^2 + p02*y^2',...
%     'dependent',{'z'},'independent',{'x', 'y'},...
%     'coefficients',{'p00','p10','p01','p20', 'p02'});
% mdl = fit([y; z]', y'-x',ft)
% 
% ft = fittype('p00 + p10*x + p01*y + p11*x*y + p02*y^2',...
%     'dependent',{'z'},'independent',{'x', 'y'},...
%     'coefficients',{'p00','p10','p01', 'p11', 'p02'});
% mdl = fit([y; z]', y'-x',ft)



% ft = fittype('a - b^x * (c1*y + c2*y^2 + c3*y^3)',...
%     'dependent',{'z'},'independent',{'x', 'y'},...
%     'coefficients',{'a', 'b', 'c1', 'c2', 'c3'});
% options = fitoptions('Method', 'NonlinearLeastSquares', ...
%                      'Lower',[-Inf 0 -Inf  -Inf], 'Upper',[Inf 0.5 Inf  Inf]);
% mdl = fit([y; z]', y'-x',ft, options);

% ft = fittype('a - b^x * (c1*y + c2*y^2 )',...
%     'dependent',{'z'},'independent',{'x', 'y'},...
%     'coefficients',{'a', 'b', 'c1', 'c2'});
% options = fitoptions('Method', 'NonlinearLeastSquares', ...
%                      'Lower',[-Inf 0.3 -Inf  0.001], 'Upper',[Inf 0.5 Inf  Inf]);
% mdl = fit([y; z]', y'-x',ft, options);

% 
% ft = fittype('a - b1^x + c1*(y - c2).^2 + c3*y',...
%     'dependent',{'z'},'independent',{'x', 'y'},...
%     'coefficients',{'a', 'b1', 'c1', 'c2', 'c3'});
% 
% ft = fittype('a - b1^x * ( c1*(y - c2).^2 + c3*y)',...
%     'dependent',{'z'},'independent',{'x', 'y'},...
%     'coefficients',{'a', 'b1', 'c1', 'c2', 'c3'});
% 
% % ft = fittype('a*(x - b).^2 + c*x + d', 'dependent', 'z', 'independent', 'x', 'coefficients', {'a', 'b', 'c', 'd'});
% % ft = fittype('-a^x + b', 'dependent', 'z', 'independent', 'x', 'coefficients', {'a', 'b'});
% options = fitoptions('Method', 'NonlinearLeastSquares', ...
%                          'Lower',[-Inf 0 -Inf 10 -Inf], 'Upper',[Inf 0.6 Inf 60 Inf]);
% mdl = fit([y; z]', y'-x',ft, options)




% ft = fittype('a - b1^x * ( c1*(y - c2)^2 + c3*y)',...
%     'dependent',{'z'},'independent',{'x', 'y'},...
%     'coefficients',{'a', 'b1', 'c1', 'c2', 'c3'});
% % ft = fittype('a*(x - b).^2 + c*x + d', 'dependent', 'z', 'independent', 'x', 'coefficients', {'a', 'b', 'c', 'd'});
% % ft = fittype('-a^x + b', 'dependent', 'z', 'independent', 'x', 'coefficients', {'a', 'b'});
% options = fitoptions('Method', 'NonlinearLeastSquares', ...
%                          'Lower',[-Inf 0.15 -0.0005 10 -Inf], 'Upper',[Inf 0.6 Inf 60 Inf]);
% mdl = fit([y; z]', y'-x',ft, options)



% ft = fittype('a - b1^x  + ( (1/x)*c1*(y - c2)^2 )',...
%     'dependent',{'z'},'independent',{'x', 'y'},...
%     'coefficients',{'a', 'b1', 'c1', 'c2'});
% options = fitoptions('Method', 'NonlinearLeastSquares', ...
%                          'Lower',[-Inf 0 0 10], 'Upper',[Inf 0.5 Inf 60 ]);
% mdl = fit([y; z]', y'-x',ft, options)

% ft = fittype('(b1*x + b2)*((y-c1)^2 + c2) + a',...
%     'dependent',{'z'},'independent',{'x', 'y'},...
%     'coefficients',{'a', 'b1', 'b2', 'c1', 'c2'});
% options = fitoptions('Method', 'NonlinearLeastSquares', ...
%                          'Lower',[-Inf -Inf 0 10 -Inf], 'Upper',[Inf 0 Inf 60 Inf]);
% ft = fittype('(b1*exp(-b2*x^2))*((y-c1)^2 + c2) + a',...
%     'dependent',{'z'},'independent',{'x', 'y'},...
%     'coefficients',{'a', 'b1', 'b2', 'c1', 'c2'});
% options = fitoptions('Method', 'NonlinearLeastSquares', ...
%                          'Lower',[-Inf 0 0 10 -Inf], 'Upper',[Inf Inf Inf 60 Inf]);
% mdl = fit([y; z]', y'-x',ft, options)



X = xtrk;
Y = swh; 
Z = swh - hs; Z = -abs(Z);
nn = find(~isnan(X) & ~isnan(Y) & ~isnan(Z)); 
X = X(nn); Y = Y(nn); Z = Z(nn);

% fun = @(a) ( (a(1).*exp(-0.5*(Y./a(2)).^2)).*((X-a(3)).^2 - a(4)) ) -Z;
% ainit = [0.0008 1 40 40^2]; alb = [-Inf 0 0 -Inf]; aub = [Inf Inf 200 Inf];
% [a] = lsqnonlin(fun,ainit, alb, aub)
% mdl = @(swh, xtrk) ( (a(1).*exp(-0.5*(swh./a(2)).^2)).*((xtrk-a(3)).^2 - a(4)) ); 
% 
% 
% 
% 
% 
% fun = @(a) ( ...
% (a(1).*exp(-0.5*(Y./a(2)).^2)).*...
% (a(3).*(exp(-0.5*((X)./a(4)).^2) - 1) ) ...
% ) -Z;
% ainit = [0.0008 1 0.0006 10]; alb = [-Inf 0 -Inf 0 ]; aub = [Inf Inf Inf 60];
% [a] = lsqnonlin(fun,ainit, alb, aub);
% mdl = @(swh, xtrk) ( (a(1).*exp(-0.5*(swh./a(2)).^2)).*(a(3).*(exp(-0.5*((xtrk)./a(4)).^2) - 1) ) ); 
% 
% 
% 
% fun = @(a) ( ...
%     (a(1).*exp(-0.5*(Y./a(2)).^2)).*...
%     (a(3).*(exp(-0.5*((X)./a(4)).^2) - 1) ) ...
%     ) -Z;
% ainit = [0.0008 1 0.0006 10]; alb = [-Inf 0 -Inf 0 ]; aub = [Inf Inf Inf 60];
% [a] = lsqnonlin(fun,ainit, alb, aub)
% mdl = @(swh, xtrk) ( (a(1).*exp(-0.5*(swh./a(2)).^2)).*(a(3).*(exp(-0.5*((xtrk)./a(4)).^2) - 1) ) ); 
% 
% % EXPONENTIAL AND GAUSSIAN
% fun = @(a) ( ...
%     (a(1).*Y.^a(2)).*...
%     (a(3).*(exp(-0.5*((X)./a(4)).^2) - 1) ) ...
%     ) -Z;
% ainit = [2 -0.9 0.0006 10]; alb = [-Inf -Inf -Inf 0 ]; aub = [Inf 0 Inf 60];
% [a] = lsqnonlin(fun,ainit, alb, aub)
% mdl = @(swh, xtrk) ( (a(1).*swh.^a(2)).*(a(3).*(exp(-0.5*((xtrk)./a(4)).^2) - 1) ) ); 
% 


% EXPONENTIAL AND PARABOLIC
fun = @(a) ( (a(1).*Y.^a(2)).*((X-a(3)).^2 - a(4)) + a(5)) -Z;
ainit = [2 -0.9 40 40^2 0]; alb = [-Inf -Inf  0 -Inf -Inf]; aub = [Inf 0 200 Inf Inf];
[a] = lsqnonlin(fun,ainit, alb, aub, opts);
mdl = @(swh, xtrk) ( (a(1).*swh.^a(2)).*((xtrk-a(3)).^2 - a(4)) +  a(5)); 



end