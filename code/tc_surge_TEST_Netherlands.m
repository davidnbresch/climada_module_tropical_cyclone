function hazard=tc_surge_TEST_Netherlands
% climada
% NAME:
%   tc_surge_TEST
% PURPOSE:
%   See tc_surge_TEST for the real TEST environment
%
%   A kind of silly test, just to see whether it all works with just one or
%   just ten tracks...
%
%   TEST the tropical cyclone (TC) storm surge (TS) raw hazard creation
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
% CALLING SEQUENCE:
%   tc_surge_TEST(force_recalc_ts)
% EXAMPLE:
%   tc_surge_TEST
% INPUTS:
% OPTIONAL INPUT PARAMETERS:
% OUTPUTS:
%   writes a couple files, such as entity, assets, bathymetry and a
%       hazard event set
% MODIFICATION HISTORY:
% David N. Bresch, david.bresch@gmail.com, 20140421
% David N. Bresch, david.bresch@gmail.com, 20150819, climada_global.centroids_dir introduced
%-

hazard=[]; % init output
hazard_tc=[]; % init output

global climada_global
if ~climada_init_vars,return;end % init/import global variables

% PARAMETERS
%
module_data_dir=[fileparts(fileparts(mfilename('fullpath'))) filesep 'data'];
%
% in essence, only the TEST country, TEST location, TEST_probabilistic
% and the TC track set needs to be defined, 
% all further parameters below should be working
TEST_country_name='Netherlands';
%
TEST_location=''; % set TEST_location='' to omit labeling
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
centroids_file=[climada_global.centroids_dir filesep TEST_country_name '_centroids.mat'];
% if the centroids are generated in the present code, the entitity is also stored (not needed for this TEST)
entity_file=   [climada_global.data_dir filesep 'entities' filesep TEST_country_name '_entity.mat'];
%
% 2) tropical cyclone (TC) track hard-wired below
% -----------------------------------------------
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

% 1) read the centroids
% ---------------------
if exist(centroids_file,'file')
    load(centroids_file) % load centroids
else
    % invoke the GDP_entity moduke to generate centroids and entity
    TEST_country_name_tmp=TEST_country_name;
    if strcmp(TEST_country_name,'Vietnam'),TEST_country_name_tmp='Viet Nam';end
    [centroids,entity]=climada_create_GDP_entity(TEST_country_name_tmp);
    save(centroids_file,'centroids');
    % note: assets are not needed for the TEST, but since it's convenient to store them for later use
    save(entity_file,'entity');
    % visualize assets on map
    climada_plot_entity_assets(entity,centroids,TEST_country_name);
end

% prep the region we need
centroids_rect=[min(centroids.lon) max(centroids.lon) min(centroids.lat) max(centroids.lat)];


% 2) create TC hazard event set
% -----------------------------
% note that if not in TEST mode, one might already have a fully
% probabilistic TC hazard evetn set hence does not need to (re)create
% in order for this TEST environment to work properly and almost 
% independent of core climada,, we (re)create the TC hazard event set here

if ~exist(hazard_set_file_tc,'file')
    
    % 2) read the tc tracks
    % ---------------------
    % hard-wire one test track
    tc_track(1).MaxSustainedWindUnit='kn';
    tc_track(1).CentralPressureUnit='mb';
    tc_track(1).ID_no=1;
    tc_track(1).orig_event_flag=1;
    tc_track(1).name='TEST01';
    tc_track(1).lon=[3    3  3    3  3    4  5    6];
    tc_track(1).lat=[54.5 54 53.5 53 52.5 52 51.5 51];
    tc_track(1).MaxSustainedWind=tc_track(1).lon*0+120; % very artificial
    tc_track(1).CentralPressure=tc_track(1).lon*0+1000; % kind of dummy
    tc_track(1).yyyy=tc_track(1).lon*0+2014;
    tc_track(1).mm=tc_track(1).lon*0+1;
    tc_track(1).dd=1:length(tc_track(1).lon);
    tc_track(1).hh=12;
    tc_track(1).TimeStep=tc_track(1).lon*0+24;
    tc_track(1).datenum = datenum(tc_track(1).yyyy,tc_track(1).mm,tc_track(1).dd)+tc_track(1).hh/24;
    
    tc_track=climada_tc_equal_timestep(tc_track); % see code, make higher time resolution
    
    % plot the TEST track
    if exist(entity_file,'file')
        load(entity_file) % load entity
        climada_plot_entity_assets(entity,centroids,TEST_country_name);
        hold on;
        %plot(tc_track(1).lon,tc_track(1).lat,'-r','LineWidth',5);
        plot(tc_track(1).lon,tc_track(1).lat,'Or','MarkerSize',5);
    end
    
    if TEST_probabilistic
        tc_track=climada_tc_random_walk(tc_track); % overwrite
    end
    
    hazard = climada_tc_hazard_set(tc_track, hazard_set_file_tc, centroids);
else
    fprintf('loading TC wind hazard set from %s\n',hazard_set_file_tc);
    load(hazard_set_file_tc);
end

fprintf('TC: max(max(hazard.intensity))=%f\n',full(max(max(hazard.intensity)))); % a kind of easy check

% show biggest TC event
[~,max_tc_pos]=max(sum(hazard.intensity,2)); % the maximum TC intensity

main_fig=figure('Name','tc surge raw TEST','Position',[89 223 1014 413],'Color',[1 1 1]);
subplot(1,2,1)
values=full(hazard.intensity(max_tc_pos,:)); % get one TC footprint
centroids.lon=hazard.lon; % as the gridding routine needs centroids
centroids.lat=hazard.lat;
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
title(sprintf('windfield [m/s] (%i)',max_tc_pos));

% up to here, hazard contains the tropical cyclone (TC) hazard event set

% call the CORE code
% ==================
% hazard on input: the tropical cyclone (TC) hazard event set
% hazard on output: the storm surge (TS) hazard event set
hazard=climada_ts_hazard_set(hazard,hazard_set_file_ts);

% show biggest TC event
figure(main_fig);
subplot(1,2,2)
values=full(hazard.intensity(max_tc_pos,:)); % get one tc footprint
centroids.lon=hazard.lon; % as the gridding routine needs centroids
centroids.lat=hazard.lat;
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

return
