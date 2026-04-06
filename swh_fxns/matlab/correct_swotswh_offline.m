function [mdl] = correct_swotswh_offline(varargin)

% model syntax = mdl(y,z) = mdl(swh_karin, cross_track_distance)

%%
p = inputParser;
addParameter(p, 'data', []);
addParameter(p, 'processing', {'PGC0', 'PIC0'});
parse(p, varargin{:});
datatype = p.Results.data;
processingalgs = p.Results.processing;

%% set correction from paper parameters (utlizing both datasets)

a1 = 0.0004; a2 = -0.74; b1 = 43.39; b2 = 926.75; c = -0.08;



a = [a1 a2 b1 b2 c];

% a = [ 0.000339
%      -0.696292
%      43.532583
%     999.673946
%      -0.055497];

mdl = @(swh, xtrk) ( (a(1).*swh.^a(2)).*((xtrk-a(3)).^2 - a(4)) +  a(5)); 



end