function climada_tc_track_google_earth(tc_track, add_circles, google_earth_save)
% ----------------------------------------
%% create kml-file of all historical tc tracks 
%  lines with colors according to saffir-simpson scale
% ----------------------------------------
% NAME:
%   climada_tc_track_google_earth
% PURPOSE:
%   visualisation of historical tracks with time stamp in google
% CALLING SEQUENCE:
%   climada_tc_track_google_earth(tc_track, add_circles, google_earth_save)
% EXAMPLE:
%   climada_tc_track_google_earth
% INPUTS:
%   tc_track: track set, historical or probabilistic, but only historical
%   tracks are visualized. prompted for, if not given.
% OPTIONAL INPUT PARAMETERS:
%   add_circles:        prompted for if empty, prompt suppressed if set to 0. add
%                       circles (donuts) to every node
%   google_earth_save:  file to be saved in
%                       \climada\data\tc_tracks\tc_track....kml
% OUTPUTS:
%   kml-file, visualisation of historical tracks with time stamp in google
%   earth
% RESTRICTIONS:
%   none
% MODIFICATION HISTORY:
% Lea Mueller, muellele@gmail.com, 20110724
% Lea Mueller, muellele@gmail.com, 20150123, changed nodetime_mat to datenum
%-

global climada_global
if ~climada_init_vars, return; end
if ~exist('tc_track'         ,'var'), tc_track          = []; end
if ~exist('add_circles'      ,'var'), add_circles       = []; end
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
no_generated = length(tc_track)/no_hist;


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
 

 


%% lines with colors according to saffir-simpson scale

% waitbar 
t0       = clock;
msgstr   = sprintf('processing %i HISTORICAL tracks',length(tc_track)/no_generated);
fprintf('%s (updating waitbar with estimation of time remaining every 50th track)\n', msgstr);
h        = waitbar(0,msgstr);
% first time estimate after 10 tracks, then every 100            
mod_step = 5; 

kmlStr = [];
for i = 1:no_generated:length(tc_track)
    for node_i = 1:length(tc_track(i).lon)-1
        v       = tc_track(i).MaxSustainedWind(node_i);
        v_color = find (v < v_categories);
        v_color = v_color(1);
        kmlStr = [kmlStr ...
                  ge_plot(tc_track(i).lon(node_i:node_i+1), tc_track(i).lat(node_i:node_i+1),...
                         'name',tc_track(i).name,...
                         'lineColor',colors_(v_color,:),...
                         'timeSpanStart',datestr(tc_track(i).datenum(node_i),'yyyy-mm-ddTHH:MM:SSZ'),...
                         'timeSpanStop' ,datestr(tc_track(i).datenum(end)+6/24,'yyyy-mm-ddTHH:MM:SSZ')...
                         )];
    end
    
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


%% points (donuts) with colors according to saffir-simpson scale
%  prompt for visualisation of circles
if isempty(add_circles)
    choice = questdlg('Visualize circles for nodes?','circles');
    switch choice
    case 'Yes'
        add_circles = 1;
    case 'No'
        add_circles = 0;
    case 'Cancel'
        add_circles = 0;
    end
end

if add_circles
    % waitbar 
    t0       = clock;
    msgstr   = sprintf('processing %i HISTORICAL tracks',length(tc_track)/no_generated);
    fprintf('%s (updating waitbar with estimation of time remaining every 50th track)\n', msgstr);
    h        = waitbar(0,msgstr);
    % first time estimate after 10 tracks, then every 100            
    mod_step = 5; 
    
    fprintf('adding circles to tc tracks \n')
    
    for i = 1:no_generated:length(tc_track)
        for node_i = 1:length(tc_track(i).lon)
            v       = tc_track(i).MaxSustainedWind(node_i);
            v_color = find (v < v_categories);
            v_color = v_color(1);

            kmlStr = [kmlStr ...
                      ge_point(tc_track(i).lon(node_i), tc_track(i).lat(node_i),100,...
                      'iconURL','http://maps.google.com/mapfiles/kml/shapes/donut.png',...
                      'iconColor',colors_(v_color,:),...
                      'timeSpanStart',datestr(tc_track(i).datenum(node_i),'yyyy-mm-ddTHH:MM:SSZ'),...
                      'timeSpanStop',datestr(tc_track(i).datenum(end)+6/24,'yyyy-mm-ddTHH:MM:SSZ'))];
        end
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
end


%% prompt for google_earth_save file
if isempty(google_earth_save) % local GUI
    google_earth_save = [climada_global.data_dir filesep 'tc_tracks' filesep 'Select name to save google earth visualiation tc_tracks_counrty.kml'];
    [filename, pathname] = uiputfile(google_earth_save, 'Save google earth visualiation of tc track set set as:');
    if isequal(filename,0) || isequal(pathname,0)
        return; % cancel
    else
        google_earth_save = fullfile(pathname,filename);
    end
end

[pathstr, name, ext] = fileparts(google_earth_save);
if ~strcmp(ext,'kml')
    ext = '.kml';
end
if strcmp(pathstr,'')
    pathstr = [climada_global.data_dir filesep 'tc_tracks'];
end
google_earth_save = [fullfile(pathstr,name) ext];

fprintf('saving google earth visualisation of tc track set as\n %s\n',google_earth_save);

ge_output(google_earth_save,kmlStr)



