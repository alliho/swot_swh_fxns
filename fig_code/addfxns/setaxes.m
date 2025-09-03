function [] = setaxes(obj, varargin)
%SETAXES Sets axes for multiple subplots (axes array)
% setaxes( obj, label, dim, objtype)
if ~exist('objtype')
    objtype = 'axes';
end


if strcmp(obj(1).Type, 'figure')
%     objtype = 'figure';
    graphicsarray = get(obj, 'Children');
else
    graphicsarray = obj;
end

ai = 1;
for i=1:length(graphicsarray)
    if strcmp(graphicsarray(i).Type, objtype)
        axesarray(ai) = graphicsarray(i);
        ai = ai + 1;
    elseif strcmp(graphicsarray(i).Type, 'text')
        axesarray(ai) = graphicsarray(i);
        ai = ai + 1;
    end
    
end

% for i=1:2:numel(varargin)
%     switch lower(varargin{i})
%         case 'FontSize'
%             fontsize       = varargin{i+1};
%             for a=1:length(axesarray)
%                 set(axesarray(a), 'FontSize', fontsize)
% %                 try(axesarray(i).XLim);
% %                     xlabel(axesarray(i), label)
% %                 end 
%             end
% %         case 'ndirections'
% %             ndirections      = varargin{i+1};
% %         case 'freqround'
% %             FrequenciesRound = varargin{i+1};
%         otherwise
%             error([varargin{i} ' is not a valid property for setAllAxes2 function.']);
%     end
% end


for i=1:2:numel(varargin)
    try
        for a=1:length(axesarray)
            set(axesarray(a), varargin{i}, varargin{i+1});
        end
    catch
        try
            if strcmp(varargin{i}(2:end), 'Label')
                for a=1:length(axesarray)
                    axesarray(a).(varargin{i}).String = varargin{i+1};
                end
            elseif contains(varargin{i},'Position')
                varnm = varargin{i};
%                 vari = str2num(varnm(end-1:end-1));
                varnms = strsplit(varnm, '.');
%                 varnms{2} = varnms{2}(1:end-3);
                axesarray(a).(varnms{1}).(varnms{2}) = varargin{i+1};
            end
        catch
            error([varargin{i} ' is not a valid property for setAllAxesOpts function.']);
        end
    end
    
    
end




% 
% 
% 

% 
% 
% 
% if strcmp(dim, 'x')
%     for i=1:length(axesarray)
%         try(axesarray(i).XLim);
%             xlabel(axesarray(i), label)
%         end 
%     end
%     
% elseif strcmp(dim, 'y')
%     for i=1:length(axesarray)
%         try(axesarray(i).YLim);
%             ylabel(axesarray(i), label)
%         end 
%     end
% end
%     



end

