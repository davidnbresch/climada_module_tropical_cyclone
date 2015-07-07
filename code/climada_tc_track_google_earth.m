function climada_tc_track_google_earth(tc_track, google_earth_save)
% climada_tc_track_google_earth
% MODULE:
%   tc_hazard_advanced
% NAME:
%   climada_tc_track_google_earth
% PURPOSE:
%   visualisation of historical tracks with time stamp in google
%   create kml-file of all historical tc tracks 
%   lines with colors according to saffir-simpson scale
% CALLING SEQUENCE:
%   climada_tc_track_google_earth(tc_track, google_earth_save)
% EXAMPLE:
%   climada_tc_track_google_earth
% INPUTS:
%   tc_track: track set, historical or probabilistic, but only historical
%   tracks are visualized. prompted for, if not given.
% OPTIONAL INPUT PARAMETERS:
%   google_earth_save:  file to be saved in
%                       \climada\data\tc_tracks\tc_track....kml
% OUTPUTS:
%   kmz-file, visualisation of historical tracks with time stamp in google
%   earth
% RESTRICTIONS:
%   none
% MODIFICATION HISTORY:
% Lea Mueller, muellele@gmail.com, 20110724
% Lea Mueller, muellele@gmail.com, 20150123, changed nodetime_mat to datenum
% Lea Mueller, muellele@gmail.com, 20150130, new kml toolbox
%-

global climada_global
if ~climada_init_vars, return; end
if ~exist('tc_track'         ,'var'), tc_track          = []; end
if ~exist('google_earth_save','var'), google_earth_save = []; end


% prompt for tc_track if not given
if isempty(tc_track)
    tc_track             = [climada_global.data_dir filesep 'tc_tracks' filesep '*.mat'];
    tc_track_default     = [climada_global.data_dir filesep 'tc_tracks' filesep 'Select tc track set .mat'];
    [filename, pathname] = uigetfile(tc_track, 'Select tc track set:',tc_track_default);
    if isequal(filename,0) || isequal(pathname,0)
        return; % cancel
    else
        tc_track = fullfile(pathname,filename);
    end
end
% load the tc track set, if a filename has been passed
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


%%
if isfield(tc_track,'orig_event_flag')
    no_hist      = sum([tc_track.orig_event_flag]);
else
    no_hist = length(tc_track);
end
% comment out in case there is a special selection of tracks, and not a
% constant number of generated tracks in between the historical tracks
no_generated = length(tc_track)/no_hist;
% no_generated = 1;


% check if max sustained wind in knots
for track_i = 1:length(tc_track)
    if strcmp(tc_track(track_i).MaxSustainedWindUnit,'kn') ~= 1
        fprintf('Wind not recorded in kn, conversion to kn needed')
        return
    end
end


%saffir-simpson hurricane scale in knots
v_categories = [34 64 83 96 113 135 1000];
colors_      = [ 'ffffaa00' ; %blue
                 'ff00aa55' ; %green
                 'ff00ffff' ; %yellow
                 'ff00aaff' ; %dark yellow  
                 'ff0055ff' ; %orange  
                 'ff0000ff' ; %red
                 'ff000079' ; %dark red
                ];
 

%% prompt for google_earth_save file
if isempty(google_earth_save) % local GUI
    google_earth_save = [climada_global.data_dir filesep 'tc_tracks' filesep 'Select name to save google earth visualiation tc_tracks_counrty.kmz'];
    [filename, pathname] = uiputfile(google_earth_save, 'Save google earth visualiation of tc track set set as:');
    if isequal(filename,0) || isequal(pathname,0)
        return; % cancel
    else
        google_earth_save = fullfile(pathname,filename);
    end
end

[pathstr, name, ext] = fileparts(google_earth_save);
if ~strcmp(ext,'kmz')
    ext = '.kmz';
end
if strcmp(pathstr,'')
    pathstr = [climada_global.data_dir filesep 'tc_tracks'];
end
google_earth_save = [fullfile(pathstr,name) ext];

fprintf('saving google earth visualisation of tc track set as\n %s\n',google_earth_save); 
            


%% lines with colors according to saffir-simpson scale

% waitbar 
t0       = clock;
msgstr   = sprintf('processing %i HISTORICAL tracks',length(tc_track)/no_generated);
fprintf('%s (updating waitbar with estimation of time remaining every 50th track)\n', msgstr);
h        = waitbar(0,msgstr);
% first time estimate after 10 tracks, then every 100            
mod_step = 5; 

% google_earth_save = [climada_global.data_dir filesep 'tc_tracks' filesep 'test_North_Indian_Ocean.kml'];
k = kml(google_earth_save);
for i = 1:no_generated:length(tc_track)
    kk = k.createFolder(tc_track(i).name);
    for node_i = 1:length(tc_track(i).lon)-1
        v       = tc_track(i).MaxSustainedWind(node_i);
        v_color = find (v < v_categories);
        v_color = v_color(1);
        description_str = sprintf('%s, tc track no %i, %s - %s', tc_track(i).name, i, ...
                datestr(tc_track(i).datenum(1),'dd mmm yyyy'), datestr(tc_track(i).datenum(end),'dd mmm yyyy'));
        
        kk.plot(tc_track(i).lon(node_i:node_i+1), tc_track(i).lat(node_i:node_i+1),...
             'name',tc_track(i).name,'description',description_str,...
             'lineColor',colors_(v_color,:),...
             'timeSpanBegin',datestr(tc_track(i).datenum(node_i),'yyyy-mm-ddTHH:MM:SSZ'),...
             'timeSpanEnd' ,datestr(tc_track(i).datenum(end)+6/24,'yyyy-mm-ddTHH:MM:SSZ'));
         
        kk.point(tc_track(i).lon(node_i), tc_track(i).lat(node_i),100,...
             'description',description_str,...
             'iconURL','http://maps.google.com/mapfiles/kml/shapes/donut.png',...
             'iconScale',0.5,...
             'iconColor',colors_(v_color,:),...
             'timeSpanBegin',datestr(tc_track(i).datenum(node_i),'yyyy-mm-ddTHH:MM:SSZ'),...
             'timeSpanEnd' ,datestr(tc_track(i).datenum(end)+6/24,'yyyy-mm-ddTHH:MM:SSZ')); 
    end
    kk.point(tc_track(i).lon(node_i+1), tc_track(i).lat(node_i+1),100,...
        'description',description_str,...
        'iconURL','http://maps.google.com/mapfiles/kml/shapes/donut.png',...
        'iconScale',0.5,...
        'iconColor',colors_(v_color,:),...
        'timeSpanBegin',datestr(tc_track(i).datenum(node_i),'yyyy-mm-ddTHH:MM:SSZ'),...
        'timeSpanEnd' ,datestr(tc_track(i).datenum(end)+6/24,'yyyy-mm-ddTHH:MM:SSZ')); 
    
    if mod(i,mod_step)==0
        mod_step = 50;
        t_elapsed_track   = etime(clock,t0)/i;
        tracks_remaining  = (length(tc_track)-i)/no_generated;
        t_projected_track = t_elapsed_track*tracks_remaining;
        msgstr            = sprintf('est. %i seconds left (%i tracks)', ceil(t_projected_track), tracks_remaining);
        % update waitbar
        waitbar(i/length(tc_track),h,msgstr); 
    end
    
end

close(h); % dispose waitbar

%% open visualition in google earth
k.run





