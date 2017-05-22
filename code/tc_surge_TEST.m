%function [hazard_TS,hazard_TC]=tc_surge_TEST
% climada
% NAME:
%   tc_surge_TEST
% PURPOSE:
%   TEST the tropical cyclone (TC) storm surge (TS) hazard creation
%   1) get centroids for the test country (eg Bangladesh, see PARAMETERS)
%      if they do not exist, try to run GDP_entity in order to create them
%   2) create TC wind hazard event set
%   call climada_ts_hazard_set in order to
%      3) create bathymetry file for region
%      4) create TC surge hazard event set
%   show the result
%
%   In essence, you define the country and the code checks the generation
%   of centroids, TC and TS hazard event sets
%
%   in essence a caller for code climada_ts_hazard_set
%
%   see tc_surge_plot_3d for 3D plots of surge fields
% CALLING SEQUENCE:
%   [hazard_TS,hazard_TC]=tc_surge_TEST % if run as a function
% EXAMPLE:
%   tc_surge_TEST
% INPUTS:
% OPTIONAL INPUT PARAMETERS:
% OUTPUTS:
%   writes a couple files, such as entity, assets, bathymetry and a
%       hazard event set
%   hazard_ts: the storm surge hazard event set
%   hazard_tc: the storm wind speed hazard event set
% MODIFICATION HISTORY:
% David N. Bresch, david.bresch@gmail.com, 20140420
% David N. Bresch, david.bresch@gmail.com, 20141017 module independent of location
% David N. Bresch, david.bresch@gmail.com, 20150819, climada_global.centroids_dir introduced
% David N. Bresch, david.bresch@gmail.com, 20170522, climada_hazard_plot
%-

hazard=[]; % init output

global climada_global
if ~climada_init_vars,return;end % init/import global variables

% PARAMETERS
%
module_data_dir=[fileparts(fileparts(mfilename('fullpath'))) filesep 'data'];
%
% in essence, only the TEST country, TEST location, TEST_probabilistic
% and the TC track set needs to be defined, 
% all further parameters below should be working
TEST_country_name='Bangladesh'; tc_track_file='tracks.nio.txt';
%TEST_country_name='Mozambique'; tc_track_file='tracks.she.txt';
%
TEST_location.name='  Barisal'; % first two spaces for nicer labeling
TEST_location.longitude=90+30/60+0/3600;
TEST_location.latitude =22+48/60+0/3600;
%TEST_location=''; % set TEST_location='' to omit labeling
%
% whether we run historic only (=0, default, good enough to test) or fully probabilistic
TEST_probabilistic=1; % default=0, since fast to check
%
% see comment above, unlikely one needs to change parameters below
%
% 1) centroids for study region
% -----------------------------
% define the file with centroids (geo-locations of the points we later
% evaluate and store storm surge heights at)
% see climada_create_GDP_entity to create centroids file
% if the centroids are generated in the present code, the entitity is also stored (not needed for this TEST)
entity_file=   [climada_global.data_dir filesep 'entities' filesep TEST_country_name '_entity.mat'];
%
% 2) tropical cyclone (TC) tracks
% -------------------------------
% set UNISYS TC track data file (for info, see climada_tc_read_unisys_database)
unisys_file=   [climada_global.data_dir filesep 'tc_tracks' filesep tc_track_file];
%
% 3) bathymetry parameters set in climada_ts_hazard_set
%
% 4) surge hazard event set
% -------------------------
% define the hazard event set file to store the TEST hazard event set
hazard_set_file_tc=[climada_global.data_dir filesep 'hazards' filesep TEST_country_name '_hazard_TC.mat'];
hazard_set_file_ts=[climada_global.data_dir filesep 'hazards' filesep TEST_country_name '_hazard_TS.mat'];


% Calculations start
% ==================

entity=climada_nightlight_entity(TEST_country_name);
fprintf('> saving entity as %s\n',entity_file);
save(entity_file,'entity');
centroids.lon=entity.assets.lon;
centroids.lat=entity.assets.lat;
centroids.centroid_ID=1:length(centroids.lon);

% prep the region we need
centroids_rect=[min(centroids.lon) max(centroids.lon) min(centroids.lat) max(centroids.lat)];


% 2) create TC hazard event set
% -----------------------------
% note that if not in TEST mode, one might already have a fully
% probabilistic TC hazard evetn set hence does not need to (re)create
% in order for this TEST environment to work properly and almost
% independent of core climada,, we (re)create the TC hazard event set here

if ~exist(hazard_set_file_tc,'file')
    
    tc_track=climada_tc_read_unisys_database(unisys_file);
    
    if TEST_probabilistic
        
        if exist('climada_tc_track_wind_decay_calculate','file')
            % wind speed decay at track nodes after landfall
            [a,p_rel] = climada_tc_track_wind_decay_calculate(tc_track,1);
        else
            fprintf('NO inland decay, consider module tc_hazard_advanced\n');
        end
        
        tc_track=climada_tc_random_walk(tc_track); % overwrite
        
        if exist('climada_tc_track_wind_decay_calculate','file')
            % add the inland decay correction to all probabilistic nodes
            tc_track   = climada_tc_track_wind_decay(tc_track, p_rel,1);
        end
        
        % plot the tracks
        figure('Name','TC tracks','Color',[1 1 1]);
        hold on
        for event_i=1:length(tc_track) % plot all tracks
            plot(tc_track(event_i).lon,tc_track(event_i).lat,'-b');
        end % event_i
        % overlay historic (to make them visible, too)
         for event_i=1:length(tc_track)
            if tc_track(event_i).orig_event_flag
                plot(tc_track(event_i).lon,tc_track(event_i).lat,'-r');
            end
        end % event_i
        climada_plot_world_borders(2)
        box on
        axis equal
        axis(centroids_rect);
        xlabel('blue: probabilistic, red: historic');
        
    end
    
    % generate all the wind footprints
    hazard_TC = climada_tc_hazard_set(tc_track, hazard_set_file_tc, centroids);
    
else
    fprintf('loading TC wind hazard set from %s\n',hazard_set_file_tc);
    load(hazard_set_file_tc);
    hazard_TC=hazard;
end

fprintf('TC: max intensity=%f\n',full(max(max(hazard_TC.intensity)))); % a kind of easy check

% show biggest TC event
[~,max_tc_pos]=max(sum(hazard_TC.intensity,2)); % the maximum TC intensity

main_fig=figure('Name','tc surge TEST','Position',[89 223 1014 413],'Color',[1 1 1]);
subplot(1,2,1)
values=full(hazard_TC.intensity(max_tc_pos,:)); % get one TC footprint
centroids.lon=hazard_TC.lon; % as the gridding routine needs centroids
centroids.lat=hazard_TC.lat;
[X, Y, gridded_VALUE] = climada_gridded_VALUE(values,centroids);
contourf(X, Y, gridded_VALUE,200,'edgecolor','none')
hold on
plot(centroids.lon,centroids.lat,'.r','MarkerSize',1);
if isfield(centroids,'onLand')
    water_points=find(centroids.onLand==0);
    plot(centroids.lon(water_points),centroids.lat(water_points),'.b','MarkerSize',1);
end
box on
climada_plot_world_borders
axis equal
axis(centroids_rect);
colorbar
title(sprintf('windfield [m/s] (event %i)',max_tc_pos));
fprintf('max event %i\n',max_tc_pos);

% up to here, hazard_TC contains the tropical cyclone (TC) hazard event set

% call the CORE code
% ==================
hazard_TS=climada_ts_hazard_set(hazard_TC,hazard_set_file_ts,0,1);

fprintf('TS: max intensity=%f\n',full(max(max(hazard_TS.intensity)))); % a kind of easy check

%figure;climada_hazard_plot(hazard_TS,0);

% show biggest TS event
figure(main_fig);
subplot(1,2,2)
values=full(hazard_TS.intensity(max_tc_pos,:)); % get one tc footprint
centroids.lon=hazard_TS.lon; % as the gridding routine needs centroids
centroids.lat=hazard_TS.lat;
[X, Y, gridded_VALUE] = climada_gridded_VALUE(values,centroids);
contourf(X, Y, gridded_VALUE,200,'edgecolor','none')
hold on
plot(centroids.lon,centroids.lat,'.r','MarkerSize',1);
if isfield(centroids,'onLand')
    water_points=find(centroids.onLand==0);
    plot(centroids.lon(water_points),centroids.lat(water_points),'.b','MarkerSize',1);
end
box on
climada_plot_world_borders
axis equal
axis(centroids_rect);
colorbar
title('surgefield [m]');

if ~isempty(TEST_location)
    text(TEST_location.longitude,TEST_location.latitude,TEST_location.name)
    plot(TEST_location.longitude,TEST_location.latitude,'xk');
end

%return