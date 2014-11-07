function [centroids hazard] = climada_centroids_hazard_expand(centroids_ori, hazard_ori, tc_track)
% UNDOCUMENTED
%-

if ~exist('centroids_ori', 'var'), centroids_ori = []; end
if ~exist('hazard_ori'   , 'var'), hazard_ori    = []; end
if ~exist('tc_track'     , 'var'), tc_track      = []; end

% init global variables
global climada_global
if ~climada_init_vars,return;end

% prompt for hazard if not given
if isempty(hazard_ori) % local GUI
    hazard_ori               = [climada_global.data_dir filesep 'hazards' filesep '*.mat'];
    hazard_ori_default       = [climada_global.data_dir filesep 'hazards' filesep 'choose a hazard.mat'];
    [filename, pathname] = uigetfile(hazard_ori, 'Open existing hazard event set:',hazard_ori_default);
    if isequal(filename,0) || isequal(pathname,0)
        return; % cancel
    else
        hazard_ori = fullfile(pathname,filename);
    end
end
% load the hazard, if a filename has been passed
if ~isstruct(hazard_ori)
    hazard_ori_file = hazard_ori;
    load(hazard_ori_file)
    vars = whos('-file', hazard_ori_file);
    if ~strcmp(vars.name,'hazard_ori')
            hazard_ori = eval(vars.name);
            clear (vars.name)
    end
end

% prompt for centroids_ori if not given
if isempty(centroids_ori) % local GUI
    centroids_ori               = [climada_global.data_dir filesep 'system' filesep '*.mat'];
    centroids_ori_default       = [climada_global.data_dir filesep 'system' filesep 'choose centroids.mat'];
    [filename, pathname] = uigetfile(centroids_ori, 'Open existing centroidst:',centroids_ori_default);
    if isequal(filename,0) || isequal(pathname,0)
        return; % cancel
    else
        centroids_ori = fullfile(pathname,filename);
    end
end
% load the centroids_ori, if a filename has been passed
if ~isstruct(centroids_ori)
    centroids_ori_file = centroids_ori;
    load(centroids_ori_file)
    vars = whos('-file', centroids_ori_file);
    if ~strcmp(vars.name,'centroids_ori')
            centroids_ori = eval(vars.name);
            clear (vars.name)
    end
end

% prompt for tc_track if not given
if isempty(tc_track) % local GUI
    tc_track               = [climada_global.data_dir filesep 'tc_tracks' filesep '*.mat'];
    tc_track_default       = [climada_global.data_dir filesep 'tc_tracks' filesep 'choose track.mat'];
    [filename, pathname] = uigetfile(tc_track, 'Open existing tc track set:',tc_track_default);
    if isequal(filename,0) || isequal(pathname,0)
        return; % cancel
    else
        tc_track = fullfile(pathname,filename);
    end
end
% load the hazard, if a filename has been passed
if ~isstruct(tc_track)
    tc_track_file = tc_track;
    load(tc_track_file)
    vars = whos('-file', tc_track_file);
    if ~strcmp(vars.name,'tc_track')
            tc_track = eval(vars.name);
            clear (vars.name)
    end
end

%% CENTROIDS
% centroids_ori = centroids;
climada_figuresize(0.6,0.8);
climada_plot_world_borders
plot(centroids_ori.Longitude, centroids_ori.Latitude,'.r')
titlestr = sprintf('Original centroids, %d',size(centroids_ori.Longitude,2));
title(titlestr);


%% EXTENDED CENTROIDS
% create extended centroids, with coarse grid (1?), for nice pictures for
% wind footprint and loss footprint
% centroids    = centroids_ori;
grid_resolution  = 1; % regular grids resolution in degree
lon = [min(centroids_ori.Longitude) max(centroids_ori.Longitude)];
lat = [min(centroids_ori.Latitude ) max(centroids_ori.Latitude )];

% for US
% lon = [-100 -70];
% lat = [  20  45];

lon = lon(1):grid_resolution:lon(2);
lat = lat(1):grid_resolution:lat(2);
[lon,lat] = meshgrid(lon,lat);
no_cen                      = numel(lon);
centroids_ext.centroid_ID   = 1:no_cen;
centroids_ext.Longitude     = lon(:)';
centroids_ext.Latitude      = lat(:)';
centroids_ext.lon           = centroids_ext.Longitude;
centroids_ext.lat           = centroids_ext.Latitude;
a = cell(1,no_cen);
for a_i = 1:no_cen
    a{a_i} = 'sea';
end
centroids_ext.country_name  = a;
centroids_ext = climada_tc_on_land(centroids_ext);
centroids_ext = rmfield(centroids_ext,'lon');
centroids_ext = rmfield(centroids_ext,'lat');
centroids_ext = climada_centroids_distance_to_coast(centroids_ext, [], 1);


%% MERGE THE EXTENDED CENTROIDS WITH THE ORIGINAL CENTROIDS
centroids              = centroids_ori;
centroids.Longitude    = [centroids_ori.Longitude     centroids_ext.Longitude];
centroids.Latitude     = [centroids_ori.Latitude      centroids_ext.Latitude];
% no_cen                 = numel(lon);
centroids.centroid_ID  = 1:length(centroids.Latitude);
centroids.onLand       = [centroids_ori.onLand        centroids_ext.onLand];
centroids.dist_to_coast= [centroids_ori.dist_to_coast centroids_ext.dist_to_coast];
centroids.country_name = [centroids_ori.country_name  centroids_ext.country_name];

% visualize for check
climada_figuresize(0.6,0.8);
climada_plot_world_borders
plot(centroids.Longitude, centroids.Latitude,'.r')
hold on
plot(centroids_ori.Longitude, centroids_ori.Latitude,'.k')
titlestr = sprintf('Merge centroids, %d',size(centroids.Longitude,2));
title(titlestr);


%% EXTENDED HAZARD
% hazard_set_file = 'hazard_Japan_extended_1degree';
hazard_set_file = 'hazard_China_extended_1degree';
hazard_ext = climada_tc_hazard_set(tc_track, hazard_set_file, centroids_ext);


%% MERGE THE HAZARD WITH THE EXTENDED HAZARD
hazard = climada_hazard_merge(hazard_ori, hazard_ext);


% %% WEAKEN THE MERGED HAZARD
% hazard_wkn = climada_hazard_distance_to_coast(hazard, centroids, tc_track);



%% expand the centroids
% centroids = centroids_ori;
% grid_resolution  = 1; % regular grids resolution in degree
% lon = [-100 -70];
% lat = [  20  45];
% lon = lon(1):grid_resolution:lon(2);
% lat = lat(1):grid_resolution:lat(2);
% [lon,lat] = meshgrid(lon,lat);
% centroids.Longitude    = [centroids_ori.Longitude lon(:)'];
% centroids.Latitude     = [centroids_ori.Latitude  lat(:)'];
% no_cen                 = numel(lon);
% centroids.centroid_ID  = 1:length(centroids.Latitude);
% centroids.onLand       = [centroids_ori.onLand zeros(1,no_cen)];
% a = cell(1,no_cen);
% for a_i = 1:no_cen
%     a{a_i} = 'sea';
% end
% centroids.country_name  = [centroids_ori.country_name a];
% centroids.dist_to_coast = [centroids_ori.dist_to_coast zeros(1,no_cen)];
% 
% % centroids on land
% climada_figuresize(0.6,0.8)
% cbar = plotclr(centroids_ext.Longitude, centroids_ext.Latitude, centroids_ext.onLand,'s', 2, 0, [], [], cmap, [], 1);    



