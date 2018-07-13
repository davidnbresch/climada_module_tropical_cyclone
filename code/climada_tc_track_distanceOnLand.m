function tc_track = climada_tc_track_distanceOnLand(tc_track,time_step_size,full_recovery)

% NAME:
%   climada_tc_track_distanceOnLand
% MODULE:
%   tropical_cyclone
% PURPOSE:
%   Add field distOnLand_km to tc_track
% CALLING SEQUENCE:
%   tc_track = climada_tc_track_distanceOnLand(tc_track, full_recovery)
%
% EXAMPLE:
%   tc_track = climada_tc_track_distanceOnLand(tc_track, 1)
%
% INPUTS:
%   none, if tc_track empty prompted for
% OPTIONAL INPUT PARAMETERS:
%   time_step_size: time step size of output. 0: no change (default)
%   full_recovery: = 1, distOnLand_km is set to zero when track re-enters sea (default)
%                  = 0, distOnLand_km is reduced by distance over sea when
%                           re-entering sea
%                  = -1, distOnLand_km is perserved over sea
% OUTPUTS:
%   same structure with extra field distOnLand_km
% RESTRICTIONS:
% MODIFICATION HISTORY:
% Samuel Eberenz, eberenz@posteo.eu, 20180712, init
%-

% init global variables
global climada_global
if ~climada_init_vars, return; end

%% check inputs, and set default values
if ~exist('tc_track'      , 'var'), tc_track      = []  ; end
if ~exist('full_recovery' , 'var'), full_recovery = []  ; end
if ~exist('time_step_size', 'var'), time_step_size= []  ; end

if isempty(full_recovery), full_recovery=1; end
if isempty(time_step_size), time_step_size=0; end

% prompt for tc_track if not given
if isempty(tc_track)
    tc_track             = [climada_global.data_dir filesep 'tc_tracks' filesep '*.mat'];
    tc_track_default     = [climada_global.data_dir filesep 'tc_tracks' filesep 'Select tc track .mat'];
    [filename, pathname] = uigetfile(tc_track, 'Select tc tracks:',tc_track_default);
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
    load(tc_track_file);
end

if time_step_size
    % refine tc_tracks to time_step_size hours
    tc_track = climada_tc_equal_timestep(tc_track,time_step_size);
end
% find nodes on land and over sea
if ~isfield(tc_track,'onLand')
    tc_track = climada_tc_on_land(tc_track);
end

%% loop through all tracks
for t_i = 1:length(tc_track) % test: t_i=110;
    % init
    tc_track(t_i).distOnLand_km = tc_track(t_i).onLand*0;
    land_index_ = find(diff(tc_track(t_i).onLand) == 1)+1;
    sea_index_  = find(diff(tc_track(t_i).onLand) ==-1)+1;
    sea_index_  = [sea_index_ size(tc_track(t_i).onLand,2)];
    if min(tc_track(t_i).onLand)
        land_index_ = [1 land_index_];
    end
    
    if ~isempty(land_index_) && land_index_(1)~=sea_index_(end) % track over land before the last point?
        diff_lon = diff(tc_track(t_i).lon);
        diff_lat = diff(tc_track(t_i).lat);
        distance_km=sqrt((diff_lon.*cos(tc_track(t_i).lat(2:end)./180*pi)).^2+...
            diff_lat.^2).*111.12; % approx. conversion into km
        land_before_land = 0.5*(1+[0 tc_track(t_i).onLand(1:end-1)]).*tc_track(t_i).onLand;
        sea_before_sea = 0.5*(2-[0 tc_track(t_i).onLand(1:end-1)]).*(1-tc_track(t_i).onLand);
        %%
        for node_i=max(2,land_index_(1)):sea_index_(end)
            switch full_recovery
                case 1 % distOnLand reset to 0 when on sea
                    tc_track(t_i).distOnLand_km(node_i) = (tc_track(t_i).distOnLand_km(node_i-1)+...
                        distance_km(node_i-1)*land_before_land(node_i)) * tc_track(t_i).onLand(node_i);
                case 0 % distOnLand reduced only step by step by distance on sea
                    tc_track(t_i).distOnLand_km(node_i) = max(0,tc_track(t_i).distOnLand_km(node_i-1)+...
                        distance_km(node_i-1)*land_before_land(node_i) - distance_km(node_i-1)*sea_before_sea(node_i));
                case -1 % distOnLand perserved over sea
                    tc_track(t_i).distOnLand_km(node_i) = tc_track(t_i).distOnLand_km(node_i-1)+...
                        distance_km(node_i-1)*land_before_land(node_i);
            end
        end
    end
end
end