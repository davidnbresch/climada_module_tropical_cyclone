function climada_tr_rainrate_field_animation(tc_track, centroids,...
    aggregation,check_avi)
% plot animation of rainrate field for a specific historical or
% probabilistic storm
% NAME:
%   climada_tc_rainrate_field_animation
% PURPOSE:
%   plot animation of rainrate field for a specific historical or
%   probabilistic storm, plot is produced every aggregation time step
%   (minimum 1 hour or more)
% CALLING SEQUENCE:
%   climada_tr_rainrate_field_animation(tc_track       ,...
%                                  centroids      ,...
%                                  aggregation    ,...
%                                  check_printplot)
% EXAMPLE:
%   climada_tc_rainrate_field_animation(tc_track_prob(1226), centroids, 1, 6)
% INPUTS:
%	tc_track:           just one tc_track, tc_track_prob(1)
%   centroids:          centroid mat file
% OPTIONAL INPUT PARAMETERS:
%   aggregation:        desired timestep for plots (minimum one plot per
%   hour, can be one plot for 6 hours or more)
%   check_avi:         if set to 1 will save animation as avi-file
% OUTPUTS:
%   plot of rainrate field (footprint) for every aggreagation step
%   (minimum 1 hour) for one specific storm track
% MODIFICATION HISTORY:
% Lea Mueller, 20110603
% david.bresch@gmail.com, 20140804, GIT update
% david.bresch@gmail.com, 20170828, does not work any more
%-


global climada_global
if ~climada_init_vars,return;end % init/import global variables
if ~exist('tc_track'       ,'var'), tc_track        = []; end
if ~exist('centroids'      ,'var'), centroids       = []; end
if ~exist('aggregation'    ,'var'), aggregation     = []; end
if ~exist('check_avi'      ,'var'), check_avi       = []; end
if isempty(aggregation)           , aggregation     = 6 ; end

animation_mp4_file=...
    [climada_global.data_dir filesep 'results' filesep 'animation_movie_tr'];
video_profile='MPEG-4';

if isfield(centroids,'assets')
    % centroids are entity, copy:
    entity=centroids; clear centroids
    centroids.lon=entity.assets.lon;
    centroids.lat=entity.assets.lat;
    centroids.ID=1:length(centroids.lon);
end

%---------------------------
%% Calculations
%---------------------------
% rainrate for every hour (for every node from tc_track)
% equal timestep within this routine. silent mode on
res_one = climada_tr_rainfield(tc_track, centroids, 1, 1);
stormdate = tc_track.datenum(1);
stormname = tc_track.name;
stormname(stormname == '_') = ' ';


% aggregate wind field for specific hours (unit remains mm/s)
a            = size(res_one,1);
aggregation_count = floor(a/aggregation);

if aggregation > 1
    for i = 1:aggregation_count
        res.rainrate_aggr(i,:) = mean(res_one((i-1)*aggregation+1:i*aggregation,:));
    end
else
    res.rainrate_aggr = res_one;
end


% Prepare the new file
if check_avi
    vidObj = VideoWriter(animation_mp4_file,video_profile);
    %if ~isempty(params.frame_rate), vidObj.FrameRate = params.frame_rate; end
    open(vidObj);
end

%---------------------------
%% FIGURE
%---------------------------

fig_visible='on';
fig_handle = figure('Name','animation','visible',fig_visible,'Color',[1 1 1]);

replay    = 1;

%scale figure according to range of longitude and latitude
scale  = max(centroids.lon) - min(centroids.lon);
scale2 =(max(centroids.lon) - min(centroids.lon))/...
    (min(max(centroids.lat),60)-max(min(centroids.lat),-50));
height = 0.5;
if height*scale2 > 1.2; height = 1.2/scale2; end
fig = climada_figuresize(height,height*scale2+0.15);
set(fig,'Color',[1 1 1])

%world border and tc_track
climada_plot_world_borders(0.7)
hold on
climada_plot_tc_track_stormcategory(tc_track);
hold on
% centroids
plot(centroids.lon, centroids.lat, '+r','MarkerSize',0.8,'linewidth',0.1)
hold on

axis equal
axis([min(centroids.lon)-scale/30  max(centroids.lon)+scale/30 ...
    max(min(centroids.lat),-50)-scale/30  min(max(centroids.lat),60)+scale/30])
cmap = climada_colormap('TR');
colormap(cmap)
% gridded_max_round     = 700;
gridded_max_round     = 20;
caxis([0 gridded_max_round])

t = colorbar('YTick',[0:5:gridded_max_round]);
set(get(t,'ylabel'),'String', 'Rain rate (mm h^{-1})','fontsize',8);

xlabel('Longitude','fontsize',8)
ylabel('Latitude','fontsize',8)

set(gca,'fontsize',8)

time_ = stormdate;

while replay == 1
    gridded_max_round     = 20;
    caxis([0 gridded_max_round])
    
    for agg_i = 1:aggregation_count
        [X, Y, gridded_VALUE] = climada_gridded_VALUE(full(res.rainrate_aggr(agg_i,:)), centroids);
        %MH set vlaues lower than 0.1*unit to NaN for ploting
        gridded_VALUE(gridded_VALUE<(0.1)) = NaN;
        if sum(gridded_VALUE(:)>1)>0
            [c,h]                 = contourf(X, Y, full(gridded_VALUE),50,'edgecolor','none');
            drawnow
            climada_plot_tc_track_stormcategory(tc_track);
            plot(centroids.lon, centroids.lat, '+r','MarkerSize',0.8,'linewidth',0.1)
            time_ = stormdate + (agg_i-1)*aggregation/24;
            title([stormname ', '  datestr(time_,'dd-mmm-yy HH PM')],'fontsize',8)
            if check_avi

                currFrame   = getframe(fig_handle);
                writeVideo(vidObj,currFrame);
                
            else
                pause(0.01)
            end
            %if agg_i<aggregation_count
            delete(h)
            %end
        end
    end
    
    if check_avi
        close(vidObj);
        check_avi = [];
        fprintf('movie saved in %s\n',animation_mp4_file)
        replay=0;
    end
    
    %% relevant final foot print:
    %  rainfall sum at every centroid
    
    colormap(cmap)
    gridded_max_round     = 700;
    [X, Y, gridded_VALUE] = climada_gridded_VALUE(full(res.rainrate_aggr), centroids);
    
    %set values lower than 0.1*unit to NaN for ploting
    gridded_VALUE(gridded_VALUE<(0.1)) = NaN;
    
    [c,h]                 = contourf(X, Y, full(gridded_VALUE),[0:50:gridded_max_round],'edgecolor','none');
    climada_plot_tc_track_stormcategory(tc_track);
    plot(centroids.lon, centroids.lat, '+r','MarkerSize',0.8,'linewidth',0.1)
    title([stormname ', '  datestr(time_,'dd-mmm-yy HH PM')],'fontsize',8)
    caxis([0 gridded_max_round])
    t = colorbar('YTick',[0:100:gridded_max_round]);
    set(get(t,'ylabel'),'String', 'Rainfall sum (mm)','fontsize',8);
    set(gca,'fontsize',8)
    
%     %% ask for replay
%     if isempty(check_avi)
%         choice = questdlg('Choose your next step?','Replay and or save as animation?','replay','save as animation.avi','exit','replay');
%         switch choice
%             case 'replay'
%                 delete(h)
%                 replay    = 1;
%                 check_avi = [];
%             case 'save as animation.avi'
%                 delete(h)
%                 check_avi = 1;
%                 filename = [filesep 'results' filesep 'rainrate_animation_' stormname '_' int2str(aggregation) 'h.avi'];
%                 
%                 %mov      = avifile([climada_global.data_dir filename],'compression','none','fps',2,'quality',100);
%                 fprintf('movie saved in %s\n', filename)
%             case 'exit'
%                 close
%                 return
%         end
%     end
end % replay == 1
return

% %% edit and save colormap
% colormapeditor
% %change node pointers, etc...
% gray_blue = get(gcf,'Colormap');
% save([climada_global.data_dir '\results\mozambique\colormap_gray_blue_rate'],'gray_blue')
