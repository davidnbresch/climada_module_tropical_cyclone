function fig = climada_ED_plot(EDS, entity, Percentage_Of_Value_Flag)
% visualize Annual Expected Damage per centroid as a map
% NAME:
%   climada_ED_plot
% PURPOSE:
%   given an encoded entity (portfolio) and a hazard event set, calculate
%   the event loss set (ELS)
%   Note that the waitbar consumes quite some time, so switch it off by
%   using the climada_code_optimizer, which removes all slowing code...
% CALLING SEQUENCE:
%   ELS=climada_ELS_calc(entity,hazard,annotation_name)
% EXAMPLE:
%   ELS=climada_ELS_calc(climada_assets_encode(climada_assets_read))
% INPUTS:
%   entity: a read and encoded assets file, see climada_assets_encode(climada_assets_read)
%       > promted for if not given
%   hazard: either a hazard set (struct) or a hazard set file (.mat with a struct)
%       > promted for if not given
% OPTIONAL INPUT PARAMETERS:
%   annotation_name: a free text that will appear e.g. on plots for
%       annotation, default=''
% OUTPUTS:
%   ELS, the event loss set with:
%       reference_year: the year the losses are references to
%       event_ID(event_i): the unique ID for each event_i
%       loss(event_i): the loss amount for event_i
%       Value: the sum of allValues used in the calculation (to e.g. express
%           losses in percentage of total Value)
%       frequency(event_i): the per occurrence event frequency for each event_i
%       orig_event_flag(event_i): whether an original event (=1) or a
%           probabilistic one (=0)
%       comment: a free comment, contains time for calculation
%       hazard: itself a structure, with:
%           filename: the filename of the hazard event set
%           comment: a free comment
%       assets.filename: the filename of the assets
%       vulnerability.filename: the filename of the vulnerability
%       annotation_name: a kind of default title (sometimes empty)
% MODIFICATION HISTORY:
% David N. Bresch, david.bresch@gmail.com, 20091228
%-

global climada_global
if ~climada_init_vars,return;end % init/import global variables

% poor man's version to check arguments
if ~exist('EDS'   ,'var'), EDS    = []; end
if ~exist('entity','var'), entity = []; end
if ~exist('Percentage_Of_Value_Flag','var'),Percentage_Of_Value_Flag=0;end


% PARAMETERS
% prompt for hazard_set if not given
if isempty(entity) % local GUI
    entity=[climada_global.data_dir filesep 'entities' filesep '*.mat'];
    [filename, pathname] = uigetfile(entity, 'Select encoded entity:');
    if isequal(filename,0) || isequal(pathname,0)
        return; % cancel
    else
        entity=fullfile(pathname,filename);
    end
end
% load the hazard set, if a filename has been passed
if ~isstruct(entity)
    entity_file=entity;entity=[];
    load(entity_file);
end

if ~isfield(EDS,'ED_per_cu')
    fprintf('ED per centroid not given. Unable to proceed')
    fig = []; return
end

% prompt for hazard_set if not given
if isempty(EDS) % local GUI
    EDS=[climada_global.data_dir filesep 'results' filesep '*.mat'];
    [filename, pathname] = uigetfile(ELS, 'Select EDS for EL visualisation:');
    if isequal(filename,0) || isequal(pathname,0)
        return; % cancel
    else
        EDS=fullfile(pathname,filename);
    end
end

% load the hazard set, if a filename has been passed
if ~isstruct(EDS)
    EDS_file=EDS;EDS=[];
    load(EDS_file);
end

c_index = entity.assets.centroid_index;
if size(entity.assets.lon,2) ~= size(EDS.ED_per_cu(c_index),1)
    fprintf('Number of centroids do not match. Unable to proceed')
    fig = []; return
end

% create the figure
scale  = max(entity.assets.lon) - min(entity.assets.lon);
scale2 =(max(entity.assets.lon) - min(entity.assets.lon))/...
        (min(max(entity.assets.lat),80)-max(min(entity.assets.lat),-60));
height = 0.5;
if height*scale2 > 1.2; height = 1.2/scale2; end
ax_buffer = 3; %ax_buffer = 30;
ax_lim = [min(entity.assets.lon)-scale/ax_buffer          max(entity.assets.lon)+scale/ax_buffer ...
          max(min(entity.assets.lat),-60)-scale/ax_buffer  min(max(entity.assets.lat),80)+scale/ax_buffer];    
markersizepp = polyfit([15 62],[5 3],1);
markersize   = polyval(markersizepp,ax_lim(2) - ax_lim(1));
markersize(markersize<2) = 2;
    

fig = climada_figuresize(height,height*scale2+0.15);
if Percentage_Of_Value_Flag
    dam_TIV = EDS.ELD_per_cu(c_index) ./ entity.assets.Value' *100;
    cbar = plotclr(entity.assets.lon, entity.assets.lat, dam_TIV,'s',markersize,1,...
               [],[],[],0,1);    
    %cbar = plotclr(entity.assets.lon, entity.assets.lat, dam_TIV,'s',markersize,1,...
    %           0,max(dam_TIV),[],0,0); 
    name_str = sprintf('Expected damage \n for %s', entity.assets.filename);       
else
    miv = 1; mav = 7;
    cbar = plotclr(entity.assets.lon, entity.assets.lat, EDS.ED_per_cu(c_index),'s',markersize,1,...
               miv,mav,[],[],1);
    %cbar = plotclr(entity.assets.lon, entity.assets.lat, ELS.EL_per_cu(c_index),'s',markersize,1,...
    %           [],[],[],[],1); 
    name_str = sprintf('Expected damage \n for %s', entity.assets.filename);           
end
set(fig,'Name',name_str)  


% 
% for tc_track

% figure;plot(sort(dam_TIV),'.-')

% colormap(flipud(hot))
% cbar = plotclr(entity.assets.lon, entity.assets.lat, ELS.EL_per_cu,'s',markersize,1,...
%                [],[],[],[],1); 
%  plotclr(x,y,v, marker, markersize, colorbar_on, 
% miv, mav, map, zero_off, v_exp)          
if Percentage_Of_Value_Flag
    set(get(cbar,'ylabel'),'String', 'Expected Damage (% of TIV)' ,'fontsize',12);
else
    set(get(cbar,'ylabel'),'String', 'Expected Damage (USD) (exponential)' ,'fontsize',12);
end
hold on
box on
climada_plot_world_borders(0.5)
axis(ax_lim)
axis equal
axis(ax_lim)
title(name_str)
 