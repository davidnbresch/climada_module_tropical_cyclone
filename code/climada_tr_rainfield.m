function res=climada_tr_rainfield(tc_track,centroids,~,~,~,~)
% TC rainfield calculation (rainsum)
% NAME:
%   climada_tr_rainfield
% PURPOSE:
%   given a single TC track (lat/lon,CentralPressure,MaxSustainedWind), calculate
%   the rain sum field at locations (=centroids)
%
%   Note: this code is optimized for speed, hence assumes that tc_track is
%   free of missing data, climada_tc_equal_timestep applied and
%   MaxSustainedWind calculated.
%
%   there is still the old SLOW version climada_tr_rainfield_slow, not
%   recommended, just kep for backward compatibility
%
%   mainly called from: see climada_tr_hazard_set
% CALLING SEQUENCE:
%   res=climada_tr_rainfield(tc_track,centroids,~,~,~,~)
% EXAMPLE:
%   tc_track=climada_tc_track_load('TEST_tracks.atl_hist');
%   tc_track=climada_tc_equal_timestep(tc_track);
%   centroids=climada_centroids_load('USFL_MiamiDadeBrowardPalmBeach');
%   res=climada_tr_rainfield(tc_track(68),centroids);
%   climada_color_plot(res.rainsum,res.lon,res.lat);
% INPUTS:
%   tc_track: a structure with the single track information (length(tc_track)!=1)
%       see e.g. climada_tc_read_unisys_tc_track
%       tc_track.Azimuth and/or tc_track.Celerity calculated, if not existing
%       but climada_tc_equal_timestep mist have been run and
%       tc_track.MaxSustainedWind must exist on input
%   centroids: a structure with the centroids information (see e.g.
%       climada_centroids_read):
%       centroids.lat: the latitude of the centroids
%       centroids.lon: the longitude of the centroids
% OPTIONAL INPUT PARAMETERS:
% OUTPUTS:
%   res.rainsum: the rain fall sum [mm per storm] at all centroids
%       the single-character variables refer to the Pioneer offering circular
%       that's why we kept these short names (so one can copy the OC for
%       documentation)
%   res.lat: the latitude of the centroids (=centroids.lat)
%   res.lon: the longitude of the centroids (=centroids.lon)
%   (res.rainrate: COULD be returned, but most cases not needed, see code)
%   (res.time: time in sec, coiuld be returned, see TIMING in code)
% RESTRICTIONS:
% MODIFICATION HISTORY:
% Lea Mueller, 20110606
% Martin Heyenn 20120503
% david.bresch@gmail.com, 20140804, GIT update
% david.bresch@gmail.com, 20141020, cleanup, inreach speedup
% David N. Bresch, david.bresch@gmail.com, 20150819, climada_global.centroids_dir introduced
% David N. Bresch, david.bresch@gmail.com, 20160529, major speedup, five times faster
%-

res=[]; % init

%global climada_global
if ~climada_init_vars, return; end

if ~exist('tc_track'      ,'var'), tc_track       = []; end
if ~exist('centroids'     ,'var'), centroids      = []; end

% PARAMETERS
%
% distance (in degree) around each node we process the rainfield
dlon=3;dlat=dlon; % default=3 (=5 until 20160529), 3 deg approx 300km, enough


if isempty(tc_track),return;end
if isempty(centroids),return;end

n_track_nodes  = length(tc_track.lat);
n_centroids    = length(centroids.lat);
cos_centroids_lat = cos(centroids.lat/180*pi); % calculate once for speedup

zero_vect   =zeros(1,n_centroids); % for speedup, see further below
res.rainsum = zero_vect;
res.lat     = centroids.lat;
res.lon     = centroids.lon;

% calculate MaxSustainedWind if only CentralPressure given
if ~isfield(tc_track,'MaxSustainedWind') && isfield(tc_track,'CentralPressure')
    tc_track.MaxSustainedWind = zeros(1,n_track_nodes); % init
end
% convert to kn
switch tc_track.MaxSustainedWindUnit
    case 'km/h'
        tc_track.MaxSustainedWind = tc_track.MaxSustainedWind/1.852;
    case 'mph'
        tc_track.MaxSustainedWind = tc_track.MaxSustainedWind/1.151;
    case 'm/s'
        tc_track.MaxSustainedWind = tc_track.MaxSustainedWind/1000*60*60/1.852;
    otherwise % already kn
end;
tc_track.MaxSustainedWindUnit='kn'; % after conversion

% loop over track nodes
%t0=clock; % TIMING
for node_i = 1:n_track_nodes
    
    % just process the centroids close enough (+/-dlon resp. dlat)
    inreach=abs(res.lon-tc_track.lon(node_i))<dlon & abs(res.lat-tc_track.lat(node_i))<dlat;
    
    if any(inreach)
        
        % note the slower version uses a morte preicise distance calculation, namely
        %fRadius_km=climada_nonspheric_distance_m(res.lon(inreach),res.lat(inreach),...
        %    tc_track.lon(node_i),tc_track.lat(node_i),centroid_count,inreach);
        
        fRadius_km=zero_vect; % init
        dd=((tc_track.lon(node_i)-res.lon(inreach)).*cos_centroids_lat(inreach)).^2+...
            (tc_track.lat(node_i)-res.lat(inreach)).^2; % in km^2
        fRadius_km(inreach) = sqrt(dd)*111.12; % now in km
        
        % calculate rain rate in mm/h
        rainrate    = climada_RCLIPER(tc_track.MaxSustainedWind(node_i),inreach,fRadius_km);
        %res.rainrate(i,:)= sparse(rainrate); % uses a LOT of time, unnecessary
        res.rainsum = res.rainsum + rainrate;  %total sum of mm per wind storm
        
    end
    
end % Loop over all TC Nodes

%res.time=etime(clock,t0); % TIMING

end % climada_tr_rainfield