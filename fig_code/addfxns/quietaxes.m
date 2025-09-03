function [] = quietaxes(has, axesname)
%QUIETAXES mutes labels and tick labels on axes defined in input
% Alli Ho 06/12/2023
if strcmp(axesname, 'x')
    setAllAxesOpts(has, 'XTickLabel', '')
    setAllAxesOpts(has, 'XLabel', '')
elseif strcmp(axesname, 'y')
    setAllAxesOpts(has, 'YTickLabel', '')
    setAllAxesOpts(has, 'YLabel', '')
elseif strcmp(axesname, 'r')
    setAllAxesOpts(has, 'RTickLabel', '')
elseif strcmp(axesname, 'theta')
    setAllAxesOpts(has, 'ThetaTickLabel', '')
else
    disp(['Not valid axes name']);
end

end

