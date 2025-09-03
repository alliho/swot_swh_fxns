function [yesno] = isin(var,lims, varargin)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
% [yesno] = isin(var,lims)
% yesno = (var > lims(1) & var < lims(2));
lims = sort(lims);
yesno = (var > lims(1) & var <= lims(2));
if ~isempty(varargin)
    if strcmp(varargin{1}, 'inclusive')
        yesno = (var >= lims(1) & var <= lims(2));
    elseif strcmp(varargin{1}, 'right')
        yesno = (var > lims(1) & var <= lims(2));
    elseif strcmp(varargin{1}, 'left')
        yesno = (var >= lims(1) & var < lims(2));
    end
end
end