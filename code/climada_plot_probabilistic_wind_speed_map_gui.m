function [track_req,h,g] = climada_plot_probabilistic_wind_speed_map_gui(tc_track, track_req)
% plot historical tc track (Longitude, Latitude) in world map according to
% saffir-simpson hurrican scale. Add plot of probabilistic generated sister
% storms. Historical tracks has black lines around markers to identify as
% original track.
% NAME:
%   climada_plot_probabilistic_wind_speed_map
% PURPOSE:
%   analyse visuallly historical tc track and its generated probabilistic
%   sister storms. Check Longitude, Latitude and wind speed category
%   (saffir-simpson hurricane scale) 
% CALLING SEQUENCE:
%   climada_plot_probabilistic_wind_speed_map(tc_track)
% EXAMPLE:
%   climada_plot_probabilistic_wind_speed_map
% INPUTS:
%   tc_track: probabilistic tc track set (random walk of wind speed, 
%   longitude and latitude), wind speed in knots, nodes every six hours, if
%   not given, prompted for
% OPTIONAL INPUT PARAMETERS:
%   track_req:  number of specific historical track to be displayed with
%   its probabilistic sister storms, prompts for input 
%   p:          to print figure
%   x:          to exit
%   41:         or any other track number. will be rounded to the nearest 
%               historical track.
%   enter:      to continue.
% OUTPUTS:
%   figure, printout of figure if requested
% RESTRICTIONS:
% MODIFICATION HISTORY:
% Lea Mueller, 20110628
%-


global climada_global
if ~climada_init_vars, return; end
if ~exist('tc_track'  , 'var'), tc_track  = []; end
if ~exist('track_req' , 'var'), track_req = []; end

% if isempty(tc_track)
%     load ([climada_global.data_dir '\tc_tracks\tc_tracks_mozambique_1978_2011_southwestindian_prob_V_4480'])
%     tc_track = tc_track_prob; clear tc_track_prob
%     fprintf('\n***tc_track_prob loaded, 4480 tracks, south western Indian Ocean*** \n')                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      
% end

% prompt for probabilistic tc_track if not given
if isempty(tc_track)
    %load ([climada_global.data_dir
    %'\tc_tracks\tc_tracks_mozambique_1978_2011_southwestindian_cleaned_6h'])
    tc_track = [climada_global.data_dir filesep 'tc_tracks' filesep '*.mat'];
    [filename, pathname] = uigetfile(tc_track, 'Select PROBABILISTIC tc track set:');
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
    vars = whos('-file', tc_track_file);
    load(tc_track_file);
    if ~strcmp(vars.name,'tc_track')
        tc_track = eval(vars.name);
        clear (vars.name)
    end
end

no_hist      = sum([tc_track.orig_event_flag]);
no_generated = length(tc_track)/no_hist;
ens_size     = no_generated-1;

%check if track_req is a historical track, round to nearest historial track
if track_req
    track_req = round((track_req-1)/(ens_size+1))*(ens_size+1)+1;
end

% longitude, latitude range
lon = [tc_track(track_req+[0:ens_size]).lon];
lat = [tc_track(track_req+[0:ens_size]).lat];

x_range = [min(lon)-5 max(lon)+5];
y_range = [min(lat)-5 max(lat)+5];

hb = get(gca,'PlotBoxAspectRatio');
hb = hb(1)/hb(2);

if hb/(diff(x_range)/diff(y_range))<1
    dif     = ( diff(x_range)/hb-diff(y_range) )/2;
    y_range = [y_range(1)-dif y_range(2)+dif];
    set(gca,'xlim',x_range,'ylim',y_range)
else
    dif     = ( diff(y_range)*hb-diff(x_range) )/2;
    x_range = [x_range(1)-dif x_range(2)+dif]; 
    set(gca,'xlim',x_range,'ylim',y_range)
end
climada_plot_world_borders(0.7,[], [], 1)

h=[];
for t_i = track_req+[0:ens_size]
    tc_track_.lon = tc_track(t_i).lon(1:6:end);
    tc_track_.lat = tc_track(t_i).lat(1:6:end);
    tc_track_.MaxSustainedWind = tc_track(t_i).MaxSustainedWind(1:6:end);
    h(:,end+1) = climada_plot_tc_track_stormcategory(tc_track_,8,[]);
end
g = plot(tc_track(track_req).lon(1:6:end),tc_track(track_req).lat(1:6:end),'ok','markersize',3,'linewidth',0.7);
% g = plot(tc_track(track_req).lon,tc_track(track_req).lat,'ok','markersize',3,'linewidth',0.7);
title(['Historical track ' int2str(track_req) ' and its ' [int2str(ens_size)] ' probabilistic sister storms'])
    
    
   
   