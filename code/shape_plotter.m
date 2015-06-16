function shape_plotter(shapes,label_att,varargin)
% shape_plotter
% MODULE:
%   climada core
% NAME:
%   shape_plotter
% PURPOSE:
%   easily plot shape structs containing multiple shapes
% CALLING SEQUENCE:
%   shape_plotter(shapes,label_att,varargin)
% EXAMPLE:
%   shape_plotter(shapes,'attribute name','linewidth',2,'color','r')
% INPUTS: 
%   shapes:     struct with shape file info must have fields X and Y
%   label_att:  the fieldname containing the string with which you wish to
%               label each shape
% OPTIONAL INPUT PARAMETERS:
%   varargin:   any property value pair compatible with Matlab's plot func
% OUTPUTS:
% MODIFICATION HISTORY:
% Gilles Stassen, gillesstassen@hotmail.com, 18052015, init
%-

global climada_global
if ~climada_init_vars; return; end

if ~exist('shapes',     'var'),  return;                end
% if ~exist('name',       'var'),  name = '';             end
if ~exist('label_att',  'var'),  label_att = '';        end

if ~isstruct(shapes)
    shapes_file = shapes; shapes = [];
    [fP,fN,fE] = fileparts(shapes_file);
    
    if strcmp(fE,'.mat')
        load(shapes_file);
    elseif strcmp(fE,'.shp')
        shapes = climada_shaperead(shapes_file,0);
    else
        cprintf([1 0 0],'ERROR: invalid filetype\n')
        return;
    end
end

if ~isfield(shapes,'X') || ~isfield(shapes,'Y')
    cprintf([1 0 0],'ERROR: shapes must have attributes X and Y\n')
    return;    
end

vararg_str = '';
if ~isempty(varargin)
    for arg_i = 1: length(varargin)
        if isnumeric(varargin{arg_i})
            vararg_str = [vararg_str ',' '[' num2str(varargin{arg_i}) ']'];
        elseif ischar(varargin{arg_i})
            vararg_str = [vararg_str ',' '''' varargin{arg_i} '''']; 
        end
    end
else
    vararg_str = ',''linewidth'',1,''color'',[81 81 81]/255';
end

eval_str = ['h = plot3(shapes(shape_i).X, shapes(shape_i).Y,Z' vararg_str ');'];
hold on
% legend('-DynamicLegend');

% plot each shape in struct
for shape_i = 1:length(shapes)
    Z = ones(size(shapes(shape_i).X)).*(1000000000000);
    eval(eval_str);
    if ~isempty(label_att) && ischar(label_att)
        if isfield(shapes,'BoundingBox')
            c_lon = mean(shapes(shape_i).BoundingBox(:,1));
            c_lat = mean(shapes(shape_i).BoundingBox(:,2));
        else
            c_lon = mean(shapes(shape_i).X);
            c_lat = mean(shapes(shape_i).Y);
        end
        t = text(c_lon,c_lat,Z(1),shapes(shape_i).(label_att));
        set(t,'HorizontalAlignment','center');
    end
end
% eval_str = ['h = plot(shapes(end).X, shapes(end).Y' vararg_str ', ''DisplayName'', ''' name ''');'];
eval(eval_str)

% if ~isempty(name)
%     if ~isempty(get(legend,'HandleVisibility'))
%         legappend(h,name);
%     else
%         legend(h,name);
%     end
% end
