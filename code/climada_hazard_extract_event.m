function event = climada_hazard_extract_event(hazard,event_i)
% climada
% MODULE:
% NAME:
%   climada_hazard_extract_event
% PURPOSE:
%   easily extract one event from a hazard set for faster testing/debugging/developing
% CALLING SEQUENCE:
%   event = climada_hazard_extract_event(hazard,event_i)
% EXAMPLE:
%   event = climada_hazard_extract_event
%   event = climada_hazard_extract_event(hazard,153)
% INPUTS:
%   hazard:     standard climada hazard event set struct
%   event_i:    event number of interest, can be a vector containing
%               multiple event numbers
% OPTIONAL INPUT PARAMETERS:
% OUTPUTS:
%   event:  struct of the same format as a hazard event set, but only
%           containing the events specified by the event numbers in event_i
% MODIFICATION HISTORY:
% Gilles Stassen, gillesstassen@hotmail.com, 20150428
%-

global climada_global

if ~exist('hazard', 'var'), hazard  = [];   end
if ~exist('event_i','var'), event_i = [];   end

% prompt for hazard event set if not given
if isempty(hazard) % local GUI
    hazard_file=[climada_global.data_dir filesep 'hazards' filesep '*.mat'];
    [filename, pathname] = uigetfile(hazard_file,...
        'Select a hazard event set:');
    if isequal(filename,0) || isequal(pathname,0)
        return; % cancel
    else
        hazard_file =fullfile(pathname,filename);
    end
    load(hazard_file);
end

if isempty(event_i)
    event_i = input(sprintf('choose an event of interest (1-%i): ',length(hazard.event_ID)));
end

if ~isnumeric(event_i) || max(event_i) > length(hazard.event_ID) || any(event_i) == 0
    cprintf([1 0 0], 'ERROR: invalid event chosen')
    return
end

% find any negative event numbers (for which n = -neg_event_i shall 
% correspond to the nth largest event
neg_event_i = event_i(event_i<0);
[~,event_ii] = sort(sum(full(hazard.intensity),2));
event_i = [event_i(event_i>0) event_ii(-neg_event_i)'];

event = hazard; %init

[n_events, n_centroids] = size(event.intensity);

if n_events == n_centroids
    cprintf([1 0 0],'ERROR: hazard set has same number of events as there are centroids \n')
    cprintf([1 0 0],'\t unable to identify event related fields \n')
    return
end

% check for irrelevant fields and delete
rm_flds   = {'arr_sort' 'arr_ori_sort' 'R' 'R_ori' 'intensity_fit_ori' 'intensity_fit' 'R_fit'};
rm_ndx    = isfield(event, rm_flds);
if any(rm_ndx),   event = rmfield(event, rm_flds(rm_ndx));      end

flds = fieldnames(event);

for fld_i = 1 : length(flds)
    if numel(event.(flds{fld_i})) == n_events
        event.(flds{fld_i}) = event.(flds{fld_i})(event_i);
    end
end

event.orig_event_count  = sum(event.orig_event_flag);
event.event_count       = length(event_i);

event.intensity = event.intensity(event_i,:);

if length(event_i) > 1
    event_number_str = '';
    for i = 1: length(event_i)
        event_number_str = [event_number_str sprintf(', %i', event_i(i))];
        event.comment = sprintf('events %s of %s',event_number_str,event.comment);
    end
else
    event.comment = sprintf('event %i of %s', event_i, event.comment);
end

event.date = datestr(now);

