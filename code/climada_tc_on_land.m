function tc_track = climada_tc_on_land(tc_track, border_mask)

% Add on land variable to tc track structure for every tc track (former
% climada_tc_track_on_land)
% NAME:
%   climada_tc_on_land
% PURPOSE:
%   add on land variable (1 on land, 0 over sea) for every tc track in tc 
%   track strucuture based on world map (raster file)
%   within:  climada_tc_track_wind_decay
% CALLING SEQUENCE:
%   tc_track = climada_tc_on_land(tc_track, border_mask);
% EXAMPLE:
%   tc_track = climada_tc_on_land;
% INPUTS:
%   none, if tc_track empty prompted for
% OPTIONAL INPUT PARAMETERS:
%   tc_track: a structure with the track information including
%   tc_track.onLand (1 for on land, 0 for over sea)
% OUTPUTS:
%   same structure now including tc track variable for every tc track
% RESTRICTIONS:
% MODIFICATION HISTORY:
% Lea Mueller, 20121203
% David N. Bresch, david.bresch@gmail.com, 20140716, minor edit to catch date-line issue
%-

% init global variables
global climada_global
if ~climada_init_vars, return; end

% check inputs, and set default values
if ~exist('tc_track'       , 'var'), tc_track      = []  ; end
if ~exist('border_mask'    , 'var'), border_mask   = []  ; end

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

if isempty(border_mask)
    % set modul data directory
    modul_data_dir = [fileparts(fileparts(mfilename('fullpath'))) filesep 'data'];
    filename = [modul_data_dir filesep 'border_mask_10km.mat'];
    A = exist(filename,'file');
    if A==2
        load(filename)
    else
        fprintf('Border mask file does not exist. Unable to proceed.\n')
        return
    end
end

msgstr   = sprintf('processing %i tracks\n',length(tc_track));
h        = waitbar(0,msgstr);
mod_step = 500;
catch_count=0; % init
for t_i = 1:length(tc_track)
    tc_track(t_i).onLand = tc_track(t_i).lon*0;
    i = round((tc_track(t_i).lon - (border_mask.lon_range(1)-border_mask.resolution_x/2)) /border_mask.resolution_x)+1;
    j = round((tc_track(t_i).lat - (border_mask.lat_range(1)-border_mask.resolution_y/2)) /border_mask.resolution_y)+1;
    i(i>size(border_mask.world_mask,2)) = size(border_mask.world_mask,2);
    j(j>size(border_mask.world_mask,1)) = size(border_mask.world_mask,1);
    for n_i = 1:length(i)
        try
            tc_track(t_i).onLand(n_i) = border_mask.world_mask(j(n_i),i(n_i));
        catch
            catch_count=catch_count+1;
            tc_track(t_i).onLand(n_i) = 0;
        end
    end 
    if mod(t_i,mod_step) == 0
        mod_step          = 500;
        msgstr            = sprintf('Add onLand variable to each tc track \n%i/%i tracks',t_i, length(tc_track));
        waitbar(t_i/length(tc_track),h,msgstr); % update waitbar
    end
end
close(h)

if catch_count>0
    fprintf('WARNING: %i times dateline issue encountered, ignored\n',catch_count);
end


% figure
% imagesc(border_mask.lon_range, border_mask.lat_range, border_mask.world_mask)
% set(gca,'ydir','normal')
% 
% figure
% imagesc(border_mask.world_mask)
% set(gca,'ydir','normal')


