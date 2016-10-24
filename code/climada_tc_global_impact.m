%function climada_tc_global_impact(basin,boundary_rect,markersize,show_plots)
% climada template
% MODULE:
%   tropical_cyclone
% NAME:
%   climada_tc_global_impact
% PURPOSE:
%   show global impact of tropical cyclones
%
%   previous call: climada_tc_read_unisys_database or climada_tc_track_combine
%   next call: diverse
% CALLING SEQUENCE:
%   climada_tc_global_impact(basin)
% EXAMPLE:
%   climada_tc_global_impact('atl_hist')
% INPUTS:
%   basin: the basin name, such as 'atl_hist' or 'wpa_prob', default='atl_hist'
%       if ='all', climada_tc_track_combine is invoked to combine all
%       tracks of all basins (starting with common first year)
% OPTIONAL INPUT PARAMETERS:
%   boundary_rect: the boundary to plot [minlon maxlon minlat maxlat]
%       default is whole globe
%   markersize: the size of the 'tiles', one might need to experiment a
%       bit, as the code tries (hard) to set a reasonabls default (based on
%       resolution)
%   show_plots: if =1, show on screen, =0 only sabe to animation file (default)
% OUTPUTS:
%   plots and animation
% MODIFICATION HISTORY:
% David N. Bresch, david.bresch@gmail.com, 20161023
%-

global climada_global
if ~climada_init_vars,return;end % init/import global variables

% poor man's version to check arguments
if ~exist('basin','var'),basin='';end
if ~exist('boundary_rect','var'),boundary_rect=[];end
if ~exist('markersize','var'),markersize=[];end
if ~exist('show_plots','var'),show_plots=0;end

% TEST
basin='all';
show_plots=0;

% PARAMETERS
%
show_ocean=0; % =1, whether we color then ocean blue
show_land=0; %=1, whether we color the land grey
cat_threshold=5; % default -999 to plot all
%
animation_mp4_file=[climada_global.results_dir filesep '_TC_IMPACT_SS5.mp4'];
make_mp4=1;
%
% the size of the figure, you might check with your screen first
figure_Position=[1 5 19201 1100];
%
% set default value for basin if not given
if isempty(basin),basin='atl_hist';end
%
% the global entity
entity_file=[climada_global.entities_dir filesep 'GLOBAL_10x10.mat'];
%
country_color=[.6 .6 .6]; % light gray land color (underneath assets)
%
d_lola=1; % degrees around min/max lat/lon of assets
%
% assets coloring (used for solid colored assets, ignored for circles)
assets_cmap = makeColorMap([.6 .7 .6], [.6 .7 .9], [0 .9 0],10);
%close all;plotclr(0:10,0:10,0:10,'s',20,0,0,10*1.05,assets_cmap,1,0); % TEST colormap
%
% damage coloring (used for solid colored assets, ignored for circles)
damage_cmap = makeColorMap([0.5 .7 0],[.9 0 0],10);
%close all;plotclr(0:10,0:10,0:10,'s',20,0,0,10*1.05,damage_cmap,1,0); % TEST colormap
%
% TEST (see TEST in code below, too)
%boundary_rect=[-100 -70 20 45];
%entity_file=[climada_global.entities_dir filesep 'USA_UnitedStates_entity.mat'];
%markersize=3;

% load assets
% -----------
if exist(entity_file,'file')
    fprintf('using %s\n',entity_file);
    load(entity_file)
else
    p.restrict_Values_to_country=0;
    p.save_entity=0;
    entity=climada_nightlight_global_entity(p);
end

% get rid of all zero points for later speedup
entity.assets=climada_subarray(entity.assets,find(entity.assets.Value>0));
entity.assets.centroid_index=1:length(entity.assets.lon);

% determine the maximum damage value
max_damage_Value=log(max(entity.assets.Value)*10);

% figure the marker size
if isempty(boundary_rect)
    boundary_rect=[min(entity.assets.lon)-d_lola max(entity.assets.lon)+d_lola ...
        min(entity.assets.lat)-d_lola max(entity.assets.lat)+d_lola];
    if boundary_rect(1)<-180,boundary_rect(1)=-180;end
    if boundary_rect(2) >180,boundary_rect(2) =180;end
end

if show_plots,fig_visible='on';else fig_visible='off';end
fig_handle=figure('Name',mfilename,'Position',figure_Position,'visible',fig_visible,'Color',[1 1 1]);

if show_ocean
    fprintf('plotting ocean ...');
    
    fill([boundary_rect(1) boundary_rect(1) boundary_rect(2) boundary_rect(2)],...
        [boundary_rect(3) boundary_rect(4) boundary_rect(4) boundary_rect(3)],...
        [0.9 0.9 .99],'LineWidth',1,'FaceColor',[0.6 0.7 1],'EdgeColor',[0.6 0.7 1]) % ocean blue
    hold on
    
end % show_ocean

if show_land
    fprintf(' land ...');
    % plot land in country_color (here instead of using climada_plot_world
    % borders to speed up a bit)
    map_shape_file=climada_global.map_border_file;
    shapes=climada_shaperead(map_shape_file);
    for shape_i = 1:length(shapes)
        if isfield(shapes(shape_i),'X_ALL')
            if ~isempty(shapes(shape_i).X_ALL)
                shapes(shape_i).X=shapes(shape_i).X_ALL;
                shapes(shape_i).Y=shapes(shape_i).Y_ALL;
            end
        end
        isnan_pos=find(isnan(shapes(shape_i).X)); % find sub-shapes
        i1=1; % init
        for isnan_pos_i=1:length(isnan_pos) % plot each sub-shape without NaNs
            i2=isnan_pos(isnan_pos_i)-1;
            fill(shapes(shape_i).X(i1:i2),shapes(shape_i).Y(i1:i2),country_color,'LineWidth',1,'EdgeColor',country_color)
            i1=i2+2;
        end % isnan_pos_i
    end % shape_i
else
    climada_plot_world_borders(1);
end % show_land

axis equal
axis off
box off
set(gca,'xlim',boundary_rect(1:2),'ylim',boundary_rect(3:4));

dlon=abs(diff(boundary_rect(1:2)));
dlat=abs(diff(boundary_rect(3:4)));

if isempty(markersize)
    % a crude way to get an appropriate markersize
    markersize=max(1,15-ceil(max(dlon,dlat)));
    fprintf('markersize = %i\n',markersize);
end

fprintf(' assets ...');

asset_Value=entity.assets.Value; % to scale

%LOCAL_colorplot(entity.assets.lon,entity.assets.lat,asset_Value,assets_cmap)

plotclr(entity.assets.lon,entity.assets.lat,asset_Value,...
    's',markersize,0,0,max(asset_Value)*1.05,assets_cmap,1,0);

hold off; drawnow
fprintf(' done\n');

hold on

if strcmpi(basin,'all');
    fprintf('loading and preparing all basins:\n');
    tc_track1=climada_tc_track_load('atl_hist');
    if isempty(tc_track1)
         climada_tc_get_unisys_databases('',1);
         tc_track1=climada_tc_track_load('atl_hist');
    end
    tc_track2=climada_tc_track_load('wpa_hist');
    tc_track=climada_tc_track_combine(tc_track1,tc_track2,-1);
    tc_track2=climada_tc_track_load('epa_hist');
    tc_track=climada_tc_track_combine(tc_track ,tc_track2,-1);
    tc_track2=climada_tc_track_load('nio_hist');
    tc_track=climada_tc_track_combine(tc_track ,tc_track2,-1);
    % southern hemisphere still has the dateline issue
    %tc_track2=climada_tc_track_load('she_hist');
    %tc_track=climada_tc_track_combine(tc_track ,tc_track2,-1);
    %info=climada_tc_track_info(tc_track,1); % check plot
else
    fprintf('loading and preparing %s\n',basin);
    tc_track=climada_tc_track_load(basin);
end
tc_track=climada_tc_stormcategory(tc_track);
tc_track=climada_tc_equal_timestep(tc_track,6); %3h is good enough
fprintf('tc track preparations done\n');

% Prepare the new file
if make_mp4
    vidObj = VideoWriter(animation_mp4_file,'MPEG-4');
    open(vidObj);
    fprintf('VideoWriter to %s started\n',animation_mp4_file);
end

t0        = clock;
mod_step  = 3; % first time estimate after 10 tracks, then every 100
format_str='%s';

n_tracks=length(tc_track);
% n_tracks=2; % TEST

fprintf('plotting %i tracks\n',n_tracks);
info=cell(1,length(tc_track)); % allocate
min_yyyy= 9999;max_yyyy=-9999;

centroids.lon=entity.assets.lon;
centroids.lat=entity.assets.lat;
centroids.centroid_ID=1:length(centroids.lon);

for track_i=1:n_tracks
    %for track_i=2:2 % TEST
    
    color_cat=min(max(tc_track(track_i).category,0),5);
    
    if color_cat>=cat_threshold
        
        min_yyyy=min(min_yyyy,min(tc_track(track_i).yyyy));
        max_yyyy=max(max_yyyy,max(tc_track(track_i).yyyy));
        
        dd=min(dlon/10,dlat/10);
        h_text=text(boundary_rect(1)+dd,boundary_rect(4)-dd,...
            sprintf('%4.4i',tc_track(track_i).yyyy(1)),'FontSize',32);
        
        %tc_track(track_i).MaxSustainedWind=tc_track(track_i).MaxSustainedWind*0+100; % TEST
        
        for step_i=3:length(tc_track(track_i).lon)
            % for step_i=30:length(tc_track(track_i).lon) % TEST
            
            LineWidth=1;%if color_cat>2,LineWidth=2;end
            h_step=plot(tc_track(track_i).lon(1:step_i),tc_track(track_i).lat(1:step_i),'Color',[color_cat/7+2/7 0 0],'LineWidth',LineWidth);
            
            sub_track=climada_subarray(tc_track(track_i),1:step_i);
            
            % get the windfield
            hazard = climada_tc_hazard_set(sub_track,'NOSAVE',centroids,1,0);
            
            damage_frame=0;
            
            if sum(sum(hazard.intensity))>0
                % calculate the damage
                EDS=climada_EDS_calc(entity,hazard,'',0,2);
                
                if EDS.ED>0
                    
                    % damage plot
                    damage_Value=EDS.ED_at_centroid;
                    nz_damage_pos=find(damage_Value>0);
                    
                    %fprintf(' max damage  %f (log=%f), max_damage_Value %f\n',max(damage_Value),log(max(damage_Value)),max_damage_Value);
                    
                    damage_Value=log(damage_Value(nz_damage_pos));
                    damage_Value=min(damage_Value,max_damage_Value);
                    
                    % max(asset_Value) etc. to plot 'higher' (it is a 3D plot...)
                    plotclr(entity.assets.lon(nz_damage_pos),entity.assets.lat(nz_damage_pos),damage_Value,...
                        's',markersize,0,0,max_damage_Value,damage_cmap,1,0);
                    
                    damage_frame=1;
                    
                end % EDS.ED>0
            end % % hazard.intensity non-zero
            
            % take a frame
            if show_plots,drawnow;end
%             if damage_frame
%                 fprintf('adding frame (step %i of track %i), damage %f \n',step_i,track_i,max(EDS.ED_at_centroid));
%             else
%                 fprintf('adding frame (step %i of track %i)\n',step_i,track_i);
%             end
            
            if make_mp4
                %currFrame   = getframe(fig_handle); % inlcudes title etc.
                currFrame   = getframe(gca); % bigger frame
                writeVideo(vidObj,currFrame);
            end
            
            delete(h_step) % delete single track
            
        end % step_i
        
        delete(h_text) % delete year in upper left corner
        
        % the progress management
        if mod(track_i,mod_step)==0
            mod_step          = 10; % TEST
            t_elapsed_track   = etime(clock,t0)/track_i;
            tracks_remaining  = n_tracks-track_i;
            t_projected_sec   = t_elapsed_track*tracks_remaining;
            msgstr = sprintf('est. %3.0f sec left (%i/%i events, year %i)',...
                t_projected_sec,track_i,n_tracks,tc_track(track_i).yyyy(1));
            fprintf(format_str,msgstr); % write progress to stdout
            format_str=[repmat('\b',1,length(msgstr)) '%s']; % back to begin of line
        end
        
    end % cat_threshold
    
end % track_i
fprintf(format_str,''); % move carriage to begin of line

if make_mp4
    
    currFrame   = getframe(gca); % make sure same as above!
    writeVideo(vidObj,currFrame); % write a few frames more
    writeVideo(vidObj,currFrame);
    writeVideo(vidObj,currFrame);
    writeVideo(vidObj,currFrame);
    
    close(vidObj);
    fprintf('\n\nmovie saved as %s\n', animation_mp4_file)
end

if ~show_plots,delete(fig_handle);end

%end % climada_tc_global_impact

% below a fest way to plot, but plotclr turned out to be fast, too.
% function LOCAL_colorplot(x,y,v,miv,mav,map)
% if isempty(miv),miv=0;end
% if isempty(mav),max(v);end
% if isempty(map),map=makeColorMap([.1 .1 .1], [.1 .9 .1], [.9 .1 .9],10);end
%
% color_steps = linspace(miv,mav,size(map,1));
% color_steps = [0 color_steps];
% for nc = 2:size(map,1)
%     iv = find(v > color_steps(nc) & v<= color_steps(nc+1));
%     plot(x(iv),y(iv),'s','color',map(nc,:),'markerfacecolor',map(nc,:),'markersize',markersize,'linewidth',0.1);
% end
% iv = find(v >= mav); % values above threshold
% plot(x(iv),y(iv),'s','color',map(end,:),'markerfacecolor',map(end,:),'markersize',markersize,'linewidth',0.1);
% end % LOCAL_colorplot