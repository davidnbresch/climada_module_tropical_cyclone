function parameters=climada_tc_global_impact(parameters,test_mode)
% climada template
% MODULE:
%   tropical_cyclone
% NAME:
%   climada_tc_global_impact
% PURPOSE:
%   show global impact of tropical cyclones, save animation
%
%   previous call: climada_tc_read_unisys_database or climada_tc_track_combine
%   next call: diverse
% CALLING SEQUENCE:
%   parameters=climada_tc_global_impact(parameters,test_mode)
% EXAMPLE:
%   parameters=climada_tc_global_impact([],-1) % obtain all default values
%   climada_tc_global_impact('',1) % test and return defaults
%   p.entity_file=[climada_global.entities_dir filesep 'USA_UnitedStates_Florida_entity.mat'];
%   p.basin='atl_hist';climada_tc_global_impact(p);
% INPUTS:
% OPTIONAL INPUT PARAMETERS:
%   parameters: a structure to pass on parameters, with fields as
%       (run parameters=climada_tc_global_impact([],-1) to obtain
%       all default values)
%       basin: the basin name, such as 'atl_hist' or 'wpa_prob', default='atl_hist'
%           if ='all', climada_tc_track_combine is invoked to combine all
%           tracks of all basins (starting with common first year)
%       boundary_rect: the boundary to plot [minlon maxlon minlat maxlat]
%           default is whole globe
%       markersize: the size of the 'tiles', one might need to experiment a
%           bit, as the code tries (hard) to set a reasonabls default (based on
%           resolution)
%       d_lola: (default=1) degrees around min/max lat/lon of assets
%       assets_cmap: assets coloring (used for solid colored assets)
%           Test with close all;plotclr(0:10,0:10,0:10,'s',20,0,0,10*1.05,assets_cmap,1,0);
%       damage_cmap: damage coloring (used for solid colored assets)
%           Test with close all;plotclr(0:10,0:10,0:10,'s',20,0,0,10*1.05,damage_cmap,1,0);
%       check_plot: if =1, show on screen, =0 only save to animation file (default)
%           Do NOT set check_plot=1 except for debugging, it slows down
%           substantially
%           Set check_plot=3 to stop after plotting assets to check whether
%           e.g. markersize is fine.
%           If check_plot is negative, the routine stops after processing
%           abs(check_plot) tracks.
%       STOP_after_n_tracks: stop after processing the images for n tracks,
%           (not already after looping over n tracks). Default =NaN
%       start_year: first year to deal with, default=-9999 to start with first
%           year tracks are available
%       default_min_TimeStep: the timestep in hours between nodes 
%           (defaul =24 for speed, but =6 would be nicer)
%       verbose: =1 verbose mode, =0 not (default)
%           In any case, progress is written to stdout (time elapsed, est.)
%   test_mode: if =-1, return all default parameters, =1 real test mode
%       default =0
% OUTPUTS:
%   plots and animation
% MODIFICATION HISTORY:
% David N. Bresch, david.bresch@gmail.com, 20161023
% David N. Bresch, david.bresch@gmail.com, 20161024, substantial speedup
% David N. Bresch, david.bresch@gmail.com, 20161024, switched to parameters
%-

global climada_global
if ~climada_init_vars,return;end % init/import global variables

% poor man's version to check arguments
if ~exist('parameters','var'),parameters=struct;end
if ~exist('test_mode','var'),test_mode=0;end

% check for some parameter fields we need
if ~isfield(parameters,'basin'),parameters.basin='';end
if ~isfield(parameters,'boundary_rect'),parameters.boundary_rect=[];end
if ~isfield(parameters,'markersize'),parameters.markersize=[];end
if ~isfield(parameters,'show_ocean'),parameters.show_ocean=[];end
if ~isfield(parameters,'show_land'),parameters.show_land=[];end
if ~isfield(parameters,'cat_threshold'),parameters.cat_threshold=[];end
if ~isfield(parameters,'animation_mp4_file'),parameters.animation_mp4_file='';end
if ~isfield(parameters,'figure_Position'),parameters.figure_Position=[];end
if ~isfield(parameters,'country_color'),parameters.country_color=[];end
if ~isfield(parameters,'entity_file'),parameters.entity_file=[];end
if ~isfield(parameters,'d_lola'),parameters.d_lola=[];end
if ~isfield(parameters,'assets_cmap'),parameters.assets_cmap=[];end
if ~isfield(parameters,'damage_cmap'),parameters.damage_cmap=[];end
if ~isfield(parameters,'STOP_after_n_tracks'),parameters.STOP_after_n_tracks=[];end
if ~isfield(parameters,'start_year'),parameters.start_year=[];end
if ~isfield(parameters,'default_min_TimeStep'),parameters.default_min_TimeStep=[];end

if ~isfield(parameters,'make_mp4'),parameters.make_mp4=[];end
if ~isfield(parameters,'check_plot'),parameters.check_plot=[];end
if ~isfield(parameters,'verbose'),parameters.verbose=[];end

% set default values (see header for details)
if isempty(parameters.basin),parameters.basin='all';end
if isempty(parameters.show_ocean),parameters.show_ocean=0;end
if isempty(parameters.show_land),parameters.show_land=0;end
if isempty(parameters.cat_threshold),parameters.cat_threshold=3;end
%if isempty(parameters.figure_Position),parameters.figure_Position=[1 5 1366 668];end % larger?
if isempty(parameters.figure_Position),parameters.figure_Position=[1 5 1366 668];end % MacBookAir
if isempty(parameters.country_color),parameters.country_color=[.6 .6 .6];end % light gray land color (underneath assets)
if isempty(parameters.d_lola),parameters.d_lola=1;end
if isempty(parameters.entity_file),parameters.entity_file=...
        [climada_global.entities_dir filesep 'GLOBAL_10x10.mat'];end
if isempty(parameters.animation_mp4_file),parameters.animation_mp4_file=...
        [climada_global.results_dir filesep '_TC_IMPACT.mp4'];end
if isempty(parameters.assets_cmap),parameters.assets_cmap=...
        makeColorMap([.6 .7 .6], [.6 .7 .9], [0 .9 0],10);end
if isempty(parameters.damage_cmap),parameters.damage_cmap=...
        makeColorMap([0.5 .7 0],[.9 0 0],10);end
if isempty(parameters.STOP_after_n_tracks),parameters.STOP_after_n_tracks=NaN;end
if isempty(parameters.start_year),parameters.start_year=-9999;end
if isempty(parameters.default_min_TimeStep),parameters.default_min_TimeStep=24;end

if isempty(parameters.make_mp4),parameters.make_mp4=1;end
if isempty(parameters.check_plot),parameters.check_plot=0;end
if isempty(parameters.verbose),parameters.verbose=0;end


% PARAMETERS
%
% Here only those one will very very likely not change, others in
% parameters structure above.
%
% treat extratropical transition, to avoid unrealistic wind and damage
% fields up north
climada_global.tc.extratropical_transition=1;

if test_mode==-1
    return % returns all (default) parameters
elseif test_mode==1
    parameters.basin='atl_hist';
    parameters.cat_threshold=5;
    parameters.make_mp4=1;
    parameters.markersize=4;
    parameters.check_plot=1;
    parameters.verbose=1;
    parameters.STOP_after_n_tracks=5;
    %parameters.default_min_TimeStep=6;
    parameters.entity_file=[climada_global.entities_dir filesep ...
        'USA_UnitedStates_Florida_entity.mat'];
end

% load assets
% -----------
if exist(parameters.entity_file,'file')
    fprintf('using %s\n',parameters.entity_file);
    load(parameters.entity_file)
else
    p.restrict_Values_to_country=0;
    p.save_entity=0;
    entity=climada_nightlight_global_entity(p);
end

% get rid of all zero points for later speedup
entity.assets=climada_subarray(entity.assets,find(entity.assets.Value>0));
entity.assets.centroid_index=1:length(entity.assets.lon); % init

% determine the maximum damage value
max_damage_Value=log(max(entity.assets.Value)*10);

% figure the marker size
if isempty(parameters.boundary_rect)
    parameters.boundary_rect=[min(entity.assets.lon)-parameters.d_lola max(entity.assets.lon)+parameters.d_lola ...
        min(entity.assets.lat)-parameters.d_lola max(entity.assets.lat)+parameters.d_lola];
    if parameters.boundary_rect(1)<-180,parameters.boundary_rect(1)=-180;end
    if parameters.boundary_rect(2)> 180,parameters.boundary_rect(2)= 180;end
    if parameters.boundary_rect(3)< -90,parameters.boundary_rect(3)= -90;end
    if parameters.boundary_rect(4)>  90,parameters.boundary_rect(4)=  90;end
    
end

if parameters.check_plot,fig_visible='on';else fig_visible='off';end
fig_handle=figure('Name',mfilename,'Position',parameters.figure_Position,'visible',fig_visible,'Color',[1 1 1]);

if parameters.show_ocean
    fprintf('plotting ocean ...');
    
    fill([parameters.boundary_rect(1) parameters.boundary_rect(1) parameters.boundary_rect(2) parameters.boundary_rect(2)],...
        [parameters.boundary_rect(3) parameters.boundary_rect(4) parameters.boundary_rect(4) parameters.boundary_rect(3)],...
        [0.9 0.9 .99],'LineWidth',1,'FaceColor',[0.6 0.7 1],'EdgeColor',[0.6 0.7 1]) % ocean blue
    hold on
    
end % parameters.show_ocean

if parameters.show_land
    fprintf(' land ...');
    % plot land in parameters.country_color (here instead of using climada_plot_world
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
            fill(shapes(shape_i).X(i1:i2),shapes(shape_i).Y(i1:i2),parameters.country_color,'LineWidth',1,'EdgeColor',parameters.country_color)
            i1=i2+2;
        end % isnan_pos_i
    end % shape_i
else
    climada_plot_world_borders(1); % just borders
end % parameters.show_land

axis equal
axis off
box off
set(gca,'xlim',parameters.boundary_rect(1:2),'ylim',parameters.boundary_rect(3:4));

dlon=abs(diff(parameters.boundary_rect(1:2)));
dlat=abs(diff(parameters.boundary_rect(3:4)));

if isempty(parameters.markersize)
    % a crude way to get an appropriate markersize
    parameters.markersize=max(1,15-ceil(max(dlon,dlat)));
    fprintf('markersize = %i\n',parameters.markersize);
end

fprintf(' assets ...');

asset_Value=entity.assets.Value; % to scale

%LOCAL_colorplot(entity.assets.lon,entity.assets.lat,asset_Value,parameters.assets_cmap)

plotclr(entity.assets.lon,entity.assets.lat,asset_Value,...
    's',parameters.markersize,0,0,max(asset_Value)*1.05,parameters.assets_cmap,1,0);

hold off;drawnow
fprintf(' done\n');
hold on

if parameters.check_plot>2 && isnan(parameters.STOP_after_n_tracks)
    fprintf('STOP: returned after plotting assets, markersize=%i\n',parameters.markersize);
    return
end

% load tracks
% -----------

if strcmpi(parameters.basin,'all');
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
    fprintf('loading and preparing %s\n',parameters.basin);
    tc_track=climada_tc_track_load(parameters.basin);
end
tc_track=climada_tc_stormcategory(tc_track);
tc_track=climada_tc_equal_timestep(tc_track,parameters.default_min_TimeStep); %3h is good enough
fprintf('tc track preparations done\n');

% prepare the animation file
% --------------------------
if parameters.make_mp4
    vidObj = VideoWriter(parameters.animation_mp4_file,'MPEG-4');
    open(vidObj);
    fprintf('VideoWriter to %s started\n',parameters.animation_mp4_file);
end

t0        = clock;
mod_step  = 1; % first time estimate after 2 tracks, then every 5
format_str='%s';

n_tracks=length(tc_track);

fprintf('looping over %i tracks\n',n_tracks);
min_yyyy= 9999;max_yyyy=-9999;

centroids.lon=entity.assets.lon;
centroids.lat=entity.assets.lat;
centroids.centroid_ID=1:length(centroids.lon);

boundary_poly_x=[parameters.boundary_rect(1) parameters.boundary_rect(1) parameters.boundary_rect(2) parameters.boundary_rect(2) parameters.boundary_rect(1)];
boundary_poly_y=[parameters.boundary_rect(3) parameters.boundary_rect(4) parameters.boundary_rect(4) parameters.boundary_rect(3) parameters.boundary_rect(3)];

n_tracks_plotted=0;

for track_i=1:n_tracks
    %for track_i=2:3 % TEST
    
    color_cat=min(max(tc_track(track_i).category,0),5);
    yyyy=tc_track(track_i).yyyy(1);
    
    if color_cat>=parameters.cat_threshold && yyyy>=parameters.start_year
        
        min_yyyy=min(min_yyyy,min(tc_track(track_i).yyyy));
        max_yyyy=max(max_yyyy,min(tc_track(track_i).yyyy));
        
        dd=min(dlon/10,dlat/10);
        h_text=text(parameters.boundary_rect(1)+dd,parameters.boundary_rect(4)-dd,...
            sprintf('%4.4i',yyyy),'FontSize',32);
        
        % check for track being visible
        in=inpolygon(tc_track(track_i).lon,tc_track(track_i).lat,boundary_poly_x,boundary_poly_y);
        
        %tc_track(track_i).MaxSustainedWind=tc_track(track_i).MaxSustainedWind*0+100; % TEST
        
        if sum(in)>0
            for step_i=3:length(tc_track(track_i).lon)
                %for step_i=30:length(tc_track(track_i).lon) % TEST
                
                sub_track=climada_subarray(tc_track(track_i),1:step_i);
                
                LineWidth=1;%if color_cat>2,LineWidth=2;end
                h_step=plot(sub_track.lon,sub_track.lat,'Color',[color_cat/7+2/7 0 0],'LineWidth',LineWidth);
                
                % get the windfield
                hazard = climada_tc_hazard_set(sub_track,'NOSAVE',centroids,1,0);
                
                if sum(sum(hazard.intensity))>0
                    
                    % calculate the damage
                    EDS=climada_EDS_calc(entity,hazard,'',0,2);
                    
                    if EDS.ED>0
                        
                        % damage plot
                        damage_Value=EDS.ED_at_centroid;
                        nz_damage_pos=find(damage_Value>0);
                        
                        if parameters.verbose,fprintf(' step %i of track %i (%i): max damage  %f (log=%f), max_damage_Value %f\n',...
                                step_i,track_i,yyyy,max(damage_Value),log(max(damage_Value)),max_damage_Value);end
                        
                        damage_Value=log(damage_Value(nz_damage_pos));
                        damage_Value=min(damage_Value,max_damage_Value);
                        
                        % plot damage
                        % -----------
                        plotclr(entity.assets.lon(nz_damage_pos),entity.assets.lat(nz_damage_pos),damage_Value,...
                            's',parameters.markersize,0,0,max_damage_Value,parameters.damage_cmap,1,0);
                        
                    else
                        if parameters.verbose,fprintf(' step %i of track %i (%i): max intensity %f\n',...
                                step_i,track_i,yyyy,full(max(hazard.intensity)));end
                    end % EDS.ED>0
                else
                    if parameters.verbose,fprintf(' step %i of track %i (%i)\n',step_i,track_i,yyyy);end
                end % % hazard.intensity non-zero
                
                % take a frame
                if parameters.check_plot,drawnow;end % really slowing down
                if parameters.make_mp4
                    %currFrame   = getframe(fig_handle); % inlcudes title etc.
                    currFrame   = getframe(gca); % bigger frame
                    % frame width and height need to be a multiple of two
                    if mod(size(currFrame.cdata,1),2),currFrame.cdata=currFrame.cdata(1:end-1,:,:);end
                    if mod(size(currFrame.cdata,2),2),currFrame.cdata=currFrame.cdata(:,1:end-1,:);end
                    writeVideo(vidObj,currFrame);
                end % parameters.make_mp4
                
                delete(h_step) % delete single track
                
            end % step_i
            n_tracks_plotted=n_tracks_plotted+1;
            
            % the progress management
            if mod(track_i,mod_step)==0
                mod_step         = 1; % TEST
                t_elapsed        = etime(clock,t0);
                t_elapsed_track  = t_elapsed/n_tracks_plotted;
                track_fraction   = n_tracks_plotted/track_i; % sepcial, since we do not plot all
                tracks_remaining = max(0,n_tracks*track_fraction-n_tracks_plotted);
                t_projected_sec  = t_elapsed_track*tracks_remaining;
                msgstr = sprintf('elapsed %3.0f sec, est. %3.0f sec left (plotted %i of %i of total %i tracks, year %i)',...
                    t_elapsed,t_projected_sec,n_tracks_plotted,track_i,n_tracks,yyyy);
                fprintf(format_str,msgstr); % write progress to stdout
                format_str=[repmat('\b',1,length(msgstr)) '%s']; % back to begin of line
            end
            
        else
            if parameters.verbose,fprintf(' skipped track %i (%i)\n',track_i,yyyy);end
        end % in
        
        delete(h_text) % delete year in upper left corner
        
    end % parameters.cat_threshold
    
    if n_tracks_plotted>parameters.STOP_after_n_tracks,break;end
    
end % track_i
fprintf(format_str,''); % move carriage to begin of line

fprintf('\n\nplotted %i out of %i tracks for % .. %i\n',n_tracks_plotted,n_tracks,min_yyyy,max_yyyy);

if parameters.make_mp4
    currFrame   = getframe(gca); % make sure same as above!
    if mod(size(currFrame.cdata,1),2),currFrame.cdata=currFrame.cdata(1:end-1,:,:);end
    if mod(size(currFrame.cdata,2),2),currFrame.cdata=currFrame.cdata(:,1:end-1,:);end
    writeVideo(vidObj,currFrame); % write a few frames more
    writeVideo(vidObj,currFrame);
    writeVideo(vidObj,currFrame);
    writeVideo(vidObj,currFrame);
    close(vidObj);
    fprintf('\n\nmovie saved as %s\n', parameters.animation_mp4_file)
end

if test_mode==1
    text(parameters.boundary_rect(1)+dd,parameters.boundary_rect(4)-dd,...
        'TEST done','FontSize',32,'Color',[0 1 0]);
end

if ~parameters.check_plot,delete(fig_handle);end

end % climada_tc_global_impact

% below a fest way to plot, but plotclr turned out to be fast, too.
% function LOCAL_colorplot(x,y,v,miv,mav,map,markersize)
% if isempty(miv),miv=0;end
% if isempty(mav),max(v);end
% if isempty(markersize),markersize=1;end
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