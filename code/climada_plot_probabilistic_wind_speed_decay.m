

function track_req = climada_plot_probabilistic_wind_speed_decay(tc_track, track_req, single)

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
%   climada_plot_probabilistic_wind_speed_decay(tc_track)
% EXAMPLE:
%   climada_plot_probabilistic_wind_speed_decay
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
if ~exist('single'    , 'var'), single    = []; end

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

%%
if ~isfield(tc_track,'onLand')
    fprintf('onLand variable not found in tc_track, please run climada_tc_on_land\n')
    return    
end

no_hist      = sum([tc_track.orig_event_flag]);
no_generated = length(tc_track)/no_hist;
ens_size     = no_generated-1;

ylabel('Wind speed (kn)')
xlabel('Time (h)')
hold on
box on

if single == 1
    title(['Track ' int2str(track_req) ', single track only'])

else
    %check if track_req is a historical track, round to nearest historial track
    if track_req
        track_req = round((track_req-1)/(ens_size+1))*(ens_size+1)+1;
    end
    title(['Historical track ' int2str(track_req) ' and its ' [int2str(ens_size)] ' probabilistic sister storms'])
    
    for t_i = track_req+[1:ens_size]
        onLand = logical(tc_track(t_i).onLand);
        if any(~onLand)
            plot(find(~onLand), tc_track(t_i).MaxSustainedWind(~onLand),'b.','markersize',3)
        end
        if any(onLand)
            plot(find(onLand), tc_track(t_i).MaxSustainedWind(onLand),'b.')
        end
        lf = find(diff(onLand) == 1);
        %fprintf('lf: ')
        %fprintf(' %d, ',lf')
        %fprintf('\n')
        if ~isempty(lf)
            plot(lf, tc_track(t_i).MaxSustainedWind(lf),'bo','markersize',5)
        end
    end
end


onLand = logical(tc_track(track_req).onLand);
if any( ~onLand)
    g(1) = plot(find(~onLand), tc_track(track_req).MaxSustainedWind(~onLand),'k.','markersize',3);
end
if any(onLand)
    g(2) = plot(find(onLand), tc_track(track_req).MaxSustainedWind(onLand),'k.');
end
lf   = find(diff(onLand) == 1);
if ~isempty(lf)
    g(3) = plot(lf, tc_track(track_req).MaxSustainedWind(lf),'or','markersize',7);
end
ylim([0 max(tc_track(track_req).MaxSustainedWind)*1.2])
xlim([0 size(onLand,2)*1.1])
if length(g) == 3
    legend(g,'wind on sea', 'wind on land', 'landfall','location','southwest')
end

       