

% figures for advanced hazard manual

%% TRACK
t_i = 1340; %katrina, aug 2005
climada_figuresize(0.6,0.8);
climada_plot_world_borders
hold on
climada_plot_tc_track_stormcategory(tc_track(t_i), [], 1);
axis equal
titlestr = sprintf('%s, %s - %s',strrep(tc_track(t_i).name,' ',''), ...
                        datestr(tc_track(t_i).nodetime_mat(1),'dd mmm'), ...
                        datestr(tc_track(t_i).nodetime_mat(end),'dd mmm yyyy'));
title(titlestr)
axis([-95 -60 20 45])

%% WIND DECAY, calculate
[tc_track_ p_rel] = climada_tc_track_wind_decay_calculate(tc_track, 1);
[tc_track_ p_rel] = climada_tc_track_wind_decay_calculate(tc_track(7:10:end), 1);

%% WIND DECAY, apply
climada_tc_track_wind_decay


%% CENTROIDS
climada_figuresize(0.6,0.8);
climada_plot_world_borders
plot(centroids.Longitude, centroids.Latitude,'.r')

climada_centroids_read

% centroids_ori = centroids;
centroids = centroids_ori;
grid_resolution  = 1; % regular grids resolution in degree
lon = [-100 -70];
lat = [  20  45];
lon = lon(1):grid_resolution:lon(2);
lat = lat(1):grid_resolution:lat(2);
[lon,lat] = meshgrid(lon,lat);
centroids.Longitude    = [centroids_ori.Longitude lon(:)'];
centroids.Latitude     = [centroids_ori.Latitude  lat(:)'];
no_cen                 = numel(lon);
centroids.centroid_ID  = 1:length(centroids.Latitude);
centroids.onLand       = [centroids_ori.onLand zeros(1,no_cen)];
a = cell(1,no_cen);
for a_i = 1:no_cen
    a{a_i} = 'sea';
end
centroids.country_name  = [centroids_ori.country_name a];
centroids.dist_to_coast = [centroids_ori.dist_to_coast zeros(1,no_cen)];

% centroids on land
climada_figuresize(0.6,0.8)
cbar = plotclr(centroids_ext.Longitude, centroids_ext.Latitude, centroids_ext.onLand,'s', 2, 0, [], [], cmap, [], 1);    



%% DISTANCE TO COAST FOR CENTROIDS
climada_coastline_read 
coastline    = [];
check_figure = 1;
centroids    = climada_centroids_distance_to_coast(centroids, coastline, check_figure);


%% WIND FOOTPRINT
climada_figuresize(0.6,0.8);
climada_plot_world_borders
plot(centroids.Longitude, centroids.Latitude,'.r')

[res centroids] = climada_tc_windfield(tc_track(t_i), centroids, 1, 0, 1);


%% WIND FOOTPRINT PER TIMESTEP
res         = climada_tc_windfield_timestep(tc_track(t_i),centroids,1); 
aggregation = 6;
check_avi   = [];
climada_tc_windfield_animation(tc_track(t_i),centroids,aggregation,check_avi);


%% HAZARD
hazard_set_file = 'test';
hazard = climada_tc_hazard_set(tc_track(1:10), hazard_set_file, centroids);


%% HAZARD WEAKENING
tc_track = climada_tc_stormcategory(tc_track);
hazard_ori = hazard;
check_figure = 1;
hazard = climada_hazard_distance_to_coast(hazard_ori, centroids, tc_track(1:10), check_figure);

hazard_wkn = climada_hazard_distance_to_coast(hazard_unwkn, centroids, tc_track, 0);
hazard_wkn = climada_hazard_distance_to_coast_USA(hazard_unwkn, centroids, tc_track, 0);


climada_hazard_distance_to_coast_australia(hazard, centroids, tc_track(1:10), check_figure)
climada_hazard_distance_to_coast_japan(hazard, centroids, tc_track(1:10), check_figure)
climada_hazard_distance_to_coast_china(hazard, centroids, tc_track(1:10), check_figure)



%% CLIMATE CHANGE HAZARD
hazard_save_name = 'test_clim';
screw = [];
hazard_clim  =  climada_hazard_clim_scen(hazard, tc_track(1:10), hazard_save_name, 2030, screw);



%% VISUALIZE FOOTPRINT FROM HAZARD
% hazard = hazard_wkn1;
% hazard = hazard_unwkn;
% hazard = hazard_wkn;
hazard = hazard_clim;
% all centroids
t_i     = 13391;
% nametag = strrep(tc_track(t_i).name,' ','');
nametag = [];
climada_figuresize(0.6,0.8);
hazard_     = hazard;
hazard_.arr = hazard.arr(t_i,:);
climada_plot_footprint(hazard_, tc_track(t_i), nametag)
% plot(centroids.Longitude, centroids.Latitude,'.k')

figure
plot(sort(hazard_.arr(:)))

% only original centroids on land
nametag = [];
climada_figuresize(0.6,0.8);
hazard_.arr = hazard.arr(t_i,1:26706);
hazard_.lon = hazard.lon(1:26706);
hazard_.lat = hazard.lat(1:26706);
climada_plot_footprint(hazard_, tc_track(t_i), nametag)

% only extended centroids, coarse grid
nametag = [];
climada_figuresize(0.6,0.8);
hazard_.arr = hazard.arr(t_i,26707:end);
hazard_.lon = hazard.lon(26707:end);
hazard_.lat = hazard.lat(26707:end);
climada_plot_footprint(hazard_, tc_track(t_i), nametag)


% distance to coast
climada_figuresize(0.6,0.8);
climada_plot_world_borders
cbar = plotclr(centroids.Longitude, centroids.Latitude, centroids.dist_to_coast, 's',4,1,10,[],[],0,0);
axis equal
plotclr(x,y,v, marker, markersize, colorbar_on, miv, mav, map, zero_off, v_exp)
    




%% EXTENDED CENTROIDS
% create extended centroids, with coarse grid (1°), for nice pictures for
% wind footprint and loss footprint
grid_resolution  = 1; % regular grids resolution in degree
lon = [-100 -70];
lat = [  20  45];
lon = lon(1):grid_resolution:lon(2);
lat = lat(1):grid_resolution:lat(2);
[lon,lat] = meshgrid(lon,lat);
no_cen                      = numel(lon);
centroids_ext.centroid_ID   = 1:no_cen;
centroids_ext.Longitude     = lon(:)';
centroids_ext.Latitude      = lat(:)';
centroids_ext.lon           = centroids_ext.Longitude;
centroids_ext.lat           = centroids_ext.Latitude;

centroids_ext = climada_tc_on_land(centroids_ext);
centroids_ext = rmfield(centroids_ext,'lon');
centroids_ext = rmfield(centroids_ext,'lat');

centroids_ext = climada_centroids_distance_to_coast(centroids_ext, [], 1);



%% EXTENDED HAZARD
hazard_set_file = 'hazard_US_extended_1degree';
hazard_ext = climada_tc_hazard_set(tc_track, hazard_set_file, centroids_ext);
hazard_ext = climada_hazard_distance_to_coast(hazard_ext, centroids_ext, tc_track);







