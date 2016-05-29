function hazard  = climada_tr_hazard_set(tc_track,hazard_set_file,centroids)
% generate hazard rain set from tc_tracks
% NAME:
%   climada_tr_hazard_set
% PURPOSE:
%   generate tropical cyclone hazard rain set
%   previous: likely climada_random_walk
%
%   CAVEAT: this is the FAST version (using parfor). For speedup reasons,
%   it does allocate the intensity array as zeros, not sparse (this way the
%   parfor loop can blindly store into) and then converts to sparse once it
%   stores it as hazard.intensity. Hence for very large computations, you
%   might run into memory problems, just search for MEMORY in the code and
%   switch to the respective lower memory usage option (still quite fast).
%
%   See climada_tr_hazard_set_slow for backward compatibility.
%
%   next: diverse, e.g. climada_EDS_calc
% CALLING SEQUENCE:
%   hazard = climada_tr_hazard_set(tc_track,hazard_set_file)
% EXAMPLE:
%   tc_track=climada_tc_track_load('TEST_tracks.atl_hist');
%   centroids=climada_centroids_load('USFL_MiamiDadeBrowardPalmBeach');
%   hazard=climada_tr_hazard_set(tc_track,'_TR_TEST_PARFOR',centroids);
% INPUTS:
% OPTIONAL INPUT PARAMETERS:
%   tc_track: a TC track structure, or a filename of a saved one
%       details: see e.g. climada_random_walk
%       > promted for if not given
%   hazard_set_file: the name of the hazard set file
%       > promted for if not given
%   centroids: the variable grid centroids (see climada_centroids_read)
%       a structure with
%           lon(1,:): the longitudes
%           lat(1,:): the latitudes
%           centroid_ID(1,:): a unique ID for each centroid, simplest: 1:length(Longitude)
%       or a file which contains the struct (saved after climada_centroids_read)
%       or an entity with entity.assets, centroids are then inferred from
%       if you select Cancel, a regular default grid is used, see hard-wired definition in code
% OUTPUTS:
%   hazard: a struct, the hazard event set, more for tests, since the
%       hazard set is stored as hazard_set_file, see code
%       lon(centroid_i): the longitude of each centroid
%       lat(centroid_i): the latitude of each centroid
%       centroid_ID(centroid_i): a unique ID for each centroid
%       peril_ID: just an ID identifying the peril, e.g. 'TC' for
%       tropical cyclone or 'ET' for extratropical cyclone
%       comment: a free comment, normally containing the time the hazard
%           event set has been generated
%       orig_years: the original years the set is based upon
%       orig_event_count: the original events
%       event_count: the total number of events in the set, including all
%           probabilistic ones, hence event_count>=orig_event_count
%       orig_event_flag(event_i): a flag for each event, whether it's an original
%           (1) or probabilistic (0) one
%       event_ID: a unique ID for each event
%       date: the creation date of the set
%       arr(event_i,centroid_i),sparse: the hazard intensity of event_i at
%           centroid_i
%       frequency(event_i): the frequency of each event
%       matrix_density: the density of the sparse array hazard.intensity
%       windfield_comment: a free comment, not in all hazard event sets
%       filename: the filename of the hazard event set (if passed as a
%           struct, this is often useful)
% MODIFICATION HISTORY:
% Lea Mueller, 20110722
% david.bresch@gmail.com, 20140804, GIT update
% David N. Bresch, david.bresch@gmail.com, 20150819, climada_global.centroids_dir introduced
% David N. Bresch, david.bresch@gmail.com, 20160529, parfor, much faster (factor 3-10)
%-

hazard = []; % init

% init global variables
global climada_global
if ~climada_init_vars,return;end

% check inputs
if ~exist('tc_track'       ,'var'), tc_track        = []; end
if ~exist('hazard_set_file','var'), hazard_set_file = []; end
if ~exist('centroids'      ,'var'), centroids       = []; end

% PARAMETERS
%
% check_plot commented out here and in climada_tc_windfield for speedup, see code

% since we store the hazard as sparse array, we need an a-priory estimation
% of it's density (only used for the MEMORY option)
hazard_arr_density    = 0.03; % 3% sparse hazard array density (estimated)

% define the reference year for this hazard set
hazard_reference_year = climada_global.present_reference_year; % default for present hazard is normally 2010


% prompt for tc_track if not given
if isempty(tc_track) % local GUI
    tc_track = [climada_global.data_dir filesep 'tc_tracks' filesep '*.mat'];
    [filename, pathname] = uigetfile(tc_track, 'Select tc track:');
    if isequal(filename,0) || isequal(pathname,0)
        return; % cancel
    else
        tc_track = fullfile(pathname,filename);
    end
end
if ~isstruct(tc_track) % load, if filename given
    tc_track_file=tc_track;tc_track=[];
    load(tc_track_file);
    vars = whos('-file', tc_track_file);
    load(tc_track_file);
    if ~strcmp(vars.name,'tc_track')
        tc_track = eval(vars.name);
        clear (vars.name)
    end
end


% prompt for hazard_set_file if not given
if isempty(hazard_set_file) % local GUI
    hazard_set_file = [climada_global.data_dir filesep 'hazards' filesep 'TRXX_hazard.mat'];
    [filename, pathname] = uiputfile(hazard_set_file, 'Save TR hazard set as:');
    if isequal(filename,0) || isequal(pathname,0)
        return; % cancel
    else
        hazard_set_file = fullfile(pathname,filename);
    end
end

% complete path, if missing
[fP,fN,fE]=fileparts(hazard_set_file);
if isempty(fP),hazard_set_file=[climada_global.data_dir filesep 'hazards' filesep fN fE];end

% prompt for centroids if not given
if isempty(centroids) % local GUI
    centroids = [climada_global.centroids_dir filesep '*.mat'];
    [filename, pathname] = uigetfile(centroids, 'Select centroids:');
    if isequal(filename,0) || isequal(pathname,0)
        % TEST centroids
        ii=0;
        for lon_i=-100:1:-50
            for lat_i=20:1:50
                ii = ii+1;
                centroids.lon(ii) = lon_i;
                centroids.lat (ii) = lat_i;
            end
        end
        centroids.centroid_ID = 1:length(centroids.lon);
        %return; % cancel
    else
        centroids = fullfile(pathname,filename);
    end
end

if isfield(centroids,'assets')
    % centroids contains in fact an entity
    entity=centroids; centroids=[]; % silly switch, but fastest
    centroids.lat =entity.assets.lat;
    centroids.lon=entity.assets.lon;
    centroids.centroid_ID=1:length(entity.assets.lon);
    % treat optional fields
    if isfield(entity.assets,'distance2coast_km'),centroids.distance2coast_km=entity.assets.distance2coast_km;end
    if isfield(entity.assets,'elevation_m'),centroids.elevation_m=entity.assets.elevation_m;end
    if isfield(entity.assets,'country_name'),centroids.country_name=entity.assets.country_name;end
    if isfield(entity.assets,'admin0_name'),centroids.admin0_name=entity.assets.admin0_name;end
    if isfield(entity.assets,'admin0_ISO3'),centroids.admin0_ISO3=entity.assets.admin0_ISO3;end
    if isfield(entity.assets,'admin1_name'),centroids.admin1_name=entity.assets.admin1_name;end
    if isfield(entity.assets,'admin1_code'),centroids.admin1_code=entity.assets.admin1_code;end
    clear entity
end

if ~isstruct(centroids) % load, if filename given
    centroids_file=centroids;centroids=[];
    % complete path, if missing
    [fP,fN,fE]=fileparts(centroids_file);
    if isempty(fP),centroids_file=[climada_global.centroids_dir filesep fN fE];end
    fprintf('centroids read from %s\n',centroids_file);
    load(centroids_file); % contains centrois as a variable
end

min_year   = tc_track  (1).yyyy  (1);
max_year   = tc_track(end).yyyy(end);
orig_years = max_year - min_year + 1;

n_centroids=length(centroids.lon);
n_tracks = length(tc_track);

hazard.lon              = centroids.lon;
hazard.lat              = centroids.lat;
hazard.centroid_ID      = centroids.centroid_ID;
hazard.peril_ID         = 'TR';
hazard.comment          = sprintf('TR torrential rain hazard event set, generated %s',datestr(now));
hazard.orig_years       = orig_years;
hazard.event_count      = n_tracks;
hazard.event_ID         = 1:n_tracks;
hazard.date             = datestr(now);
hazard.orig_event_count = 0; % init
hazard.orig_event_flag  = zeros(1,n_tracks);

% allocate the hazard array (sparse, to manage MEMORY)
% intensity = spalloc(n_tracks,n_centroids,...
%     ceil(n_tracks*n_centroids*hazard_arr_density));
intensity = zeros(n_tracks,n_centroids); % FASTER

t0       = clock;
fprintf('processing %i tracks @ %i centroids (parfor)\n',n_tracks,n_centroids);

if n_tracks>10000
    default_min_TimeStep=2; % speeds up calculation by factor 2
else
    default_min_TimeStep=climada_global.tc.default_min_TimeStep;
end
tc_track=climada_tc_equal_timestep(tc_track,default_min_TimeStep); % make equal timesteps

parfor track_i=1:n_tracks
    res                   = climada_tr_rainfield(tc_track(track_i),centroids);
    %intensity(track_i,:)  = sparse(res.rainsum); % MEMORY
    intensity(track_i,:)  = res.rainsum; % FASTER
end %track_i

hazard.intensity=sparse(intensity); % store into struct, sparse() to be safe
clear intensity % free up memory

% small lop to populate additional fields
for track_i=1:n_tracks
    hazard.orig_event_count         = hazard.orig_event_count+tc_track(track_i).orig_event_flag;
    hazard.orig_event_flag(track_i) = tc_track(track_i).orig_event_flag;
    hazard.yyyy(track_i)            = tc_track(track_i).yyyy(1);
    hazard.mm(track_i)              = tc_track(track_i).mm(1);
    hazard.dd(track_i)              = tc_track(track_i).dd(1);
    hazard.datenum(track_i)         = tc_track(track_i).datenum(1);
    hazard.name{track_i}            = tc_track(track_i).name;
end % track_i

t_elapsed = etime(clock,t0);
msgstr    = sprintf('generating %i RAIN fields took %f sec (%f sec/event)',n_tracks,t_elapsed,t_elapsed/n_tracks);
fprintf('%s\n',msgstr);

ens_size        = hazard.event_count/hazard.orig_event_count-1; % number of derived tracks per original one
event_frequency = 1/(orig_years*(ens_size+1));

hazard.frequency         = ones(1,hazard.event_count)*event_frequency; % not transposed, just regular
hazard.matrix_density    = nnz(hazard.intensity)/numel(hazard.intensity);
hazard.windfield_comment = msgstr;
hazard.filename          = hazard_set_file;
hazard.reference_year    = hazard_reference_year;

fprintf('saving TR rain hazard set as %s\n',hazard_set_file);
save(hazard_set_file,'hazard')

end % climada_tr_hazard_set