function climada_plot_lossfootprint(event_loss, centroids, tc_track, nametag)
% TC footprint figure
% NAME:
%   climada_plot_windfield
% PURPOSE:
% create footprint figure
% CALLING SEQUENCE:
%   [contr t_handle] = climada_plot_windfield(hazard, tc_track, track_no)
% EXAMPLE:
%   climada_plot_windfield
% INPUTS:
%   hazard: hazard.arr with wind intensities per centroid
%   tc_track: a structure with the track information:
%   track_no: number of track to show footprint
% OPTIONAL INPUT PARAMETERS:
% OUTPUTS:
%   figure with footprint
% RESTRICTIONS:
% MODIFICATION HISTORY:
% Lea Müller, david.bresch@gmail.com, 20121205
%-

global climada_global
if ~climada_init_vars, return; end
if ~exist('event_loss'   , 'var'), event_loss   = []; end
if ~exist('centroids'    , 'var'), centroids    = []; end
if ~exist('tc_track'     , 'var'), tc_track     = []; end
if ~exist('nametag'      , 'var'), nametag      = ''; end

% prompt for event_loss if not given
if isempty(event_loss)
    return
end

if isempty(centroids)
    centroids               = [climada_global.system_dir filesep '*.mat'];
    centroids_default       = [climada_global.system_dir filesep 'Select centroids .mat'];
    [filename, pathname] = uigetfile(centroids, 'Select hazard set:',centroids_default);
    if isequal(filename,0) || isequal(pathname,0)
        return; % cancel
    else
        centroids = fullfile(pathname,filename);
    end
end
% load the hazard set, if a filename has been passed
if ~isstruct(centroids)
    centroids_file = centroids;
    centroids      = [];
    vars = whos('-file', centroids_file);
    load(centroids_file);
    if ~strcmp(vars.name,'centroids')
        centroids = eval(vars.name);
        clear (vars.name)
    end
end
% prompt for tc_track if not given
if isempty(tc_track)
    tc_track             = [climada_global.data_dir filesep 'tc_tracks' filesep '*.mat'];
    tc_track_default     = [climada_global.data_dir filesep 'tc_tracks' filesep 'Select PROBABILISTIC tc track .mat'];
    [filename, pathname] = uigetfile(tc_track, 'Select PROBABILISTIC tc track set:',tc_track_default);
    if isequal(filename,0) || isequal(pathname,0)
        return; % cancel
    else
        tc_track = fullfile(pathname,filename);
    end
end
% load the tc track set, if a filename has been passed
if ~isstruct(tc_track)
    tc_track_file = tc_track;
    tc_track      = [];
    vars = whos('-file', tc_track_file);
    load(tc_track_file);
    if ~strcmp(vars.name,'tc_track')
        tc_track = eval(vars.name);
        clear (vars.name)
    end
    prompt   ='Type specific No. of track to print windfield [e.g. 1, 10, 34, 1011]:';
    name     =' No. of track';
    defaultanswer = {'1011'};
    answer   = inputdlg(prompt,name,1,defaultanswer);
    track_no = str2double(answer{1});
    tc_track = tc_track(track_no);
end


% set the axis
% longitude, latitude range
lon     = tc_track.lon;
lat     = tc_track.lat;
x_range = [min(lon)-5 max(lon)+5];
y_range = [min(lat)-5 max(lat)+5];
hb      = get(gca,'PlotBoxAspectRatio');
hb      = hb(1)/hb(2);
if hb/(diff(x_range)/diff(y_range))<1
    dif     = ( diff(x_range)/hb-diff(y_range) )/2;
    y_range = [y_range(1)-dif y_range(2)+dif];
    set(gca,'xlim',x_range,'ylim',y_range)
else
    dif     = ( diff(y_range)*hb-diff(x_range) )/2;
    x_range = [x_range(1)-dif x_range(2)+dif]; 
    set(gca,'xlim',x_range,'ylim',y_range)
end


    
if any(full(event_loss))
    % cmap_= [1.0000    1.0000    1.0000;
    %         0.8100    0.8100    0.8100;
    %         0.6300    0.6300    0.6300;
    %         1.0000    0.8000    0.2000;
    %         0.9420    0.6667    0.1600;
    %         0.8839    0.5333    0.1200;
    %         0.8259    0.4000    0.0800;
    %         0.7678    0.2667    0.0400;
    %         0.7098    0.1333         0];
    cmap_ = makeColorMap([255 236 139]/255, [135 206 235]/255, [39 64 139]/255,9);
    cbar  = plotclr(centroids.Longitude, centroids.Latitude, event_loss, 's',6,1,[],[],cmap_,0,1);
    set(get(cbar,'ylabel'),'String', 'Loss (USD, exponential)' ,'fontsize',8);
    
    % % create gridded values
    % [X, Y, gridded_VALUE] = climada_gridded_VALUE(full(event_loss), centroids);
    % gridded_max       = max(max(gridded_VALUE));
    % gridded_max_round = 45;
    % contourf(X, Y, full(gridded_VALUE),...
    %                 0:10:gridded_max_round,'edgecolor','none')
end
climada_plot_world_borders(0.7,[], [], 1)
climada_plot_tc_track_stormcategory(tc_track, 3, [], 0.3);
% plot(centroids.Longitude,centroids.Latitude,'.k','markersize',0.2,'linewidth',0.1)

% if any(full(windfield))
%     caxis([0 gridded_max_round])      
%     cmap_= [1.0000    1.0000    1.0000;
%             0.8100    0.8100    0.8100;
%             0.6300    0.6300    0.6300;
%             1.0000    0.8000    0.2000;
%             0.9420    0.6667    0.1600;
%             0.8839    0.5333    0.1200;
%             0.8259    0.4000    0.0800;
%             0.7678    0.2667    0.0400;
%             0.7098    0.1333         0];
%     colormap(cmap_)
%     colorbartick           = [0:10:gridded_max_round round(gridded_max)];
%     colorbarticklabel      = num2cell(colorbartick);
%     colorbarticklabel{end} = [num2str(gridded_max,'%10.2f') 'max'];
%     colorbarticklabel{end} = [int2str(gridded_max)          'max'];
%     t = colorbar('YTick',colorbartick,'yticklabel',colorbarticklabel);
%     set(get(t,'ylabel'),'String', 'Wind speed (m s^{-1})','fontsize',8);
% end
title_str = sprintf('tc track %s: %s, \t %s - %s, category %d', nametag, tc_track.name, ...
                    datestr(tc_track.nodetime_mat(1),'dd/mmm') ,datestr(tc_track.nodetime_mat(end),'dd/mmm/yyyy'),...
                    tc_track.category);               
title(title_str,'interpreter','none','fontsize',8)
% if isempty(check_printplot)
%     choice = questdlg('print?','print');
%     switch choice
%     case 'Yes'
%         check_printplot = 1;
%     case 'No'
%         check_printplot = 0;
%     case 'Cancel'
%         return
%     end
% end
% 
% if check_printplot %(>=1)   
%     foldername = [filesep 'results' filesep 'footprint_' tc_track.name '.pdf'];
%     print(fig,'-dpdf',[climada_global.data_dir foldername])
%     %close
%     fprintf('saved 1 FIGURE in folder %s \n', foldername);
% end

 
end
    
