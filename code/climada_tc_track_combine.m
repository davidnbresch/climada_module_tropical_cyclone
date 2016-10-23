function [tc_track,info]=climada_tc_track_combine(tc_track1,tc_track2,start_year)
% climada template
% MODULE:
%   tropical_cyclone
% NAME:
%   climada_tc_track_combine
% PURPOSE:
%   combine two tc_track structures, i.e merge them , such that the
%   temporal sequence is correct
%
%   previous call: climada_tc_read_unisys_database and climada_tc_track_load
%   next call: any tc code, such as climada_tc_track_info
% CALLING SEQUENCE:
%   tc_track=climada_tc_track_combine(tc_track1,tc_track2)
% EXAMPLE:
%   tc_track1=climada_tc_track_load('atl_hist');
%   tc_track2=climada_tc_track_load('wpa_hist');
%   tc_track=climada_tc_track_combine(tc_track1,tc_track2,-1);
%   info=climada_tc_track_info(tc_track,1) % check plot
%   tc_track2=climada_tc_track_load('epa_hist');
%   tc_track=climada_tc_track_combine(tc_track ,tc_track2,-1);
%   tc_track2=climada_tc_track_load('nio_hist');
%   tc_track=climada_tc_track_combine(tc_track ,tc_track2,-1);
%   tc_track2=climada_tc_track_load('she_hist');
%   tc_track=climada_tc_track_combine(tc_track ,tc_track2,-1);
% INPUTS:
%   tc_track1: a TC track strcuture, as returned from e.g. climada_tc_read_unisys_database
%   tc_track2: a second TC track strcuture
% OPTIONAL INPUT PARAMETERS:
%   start_year: a starting year, thus tracks prior to this year are ignored
%       default: empty. using all years
%       if =-1, take earliest year where both sets have tracks
% OUTPUTS:
%   tc_track: the combined tc_track structure
%   info: the (internal) info as used to compile the combined set, see code
% MODIFICATION HISTORY:
% David N. Bresch, david.bresch@gmail.com, 20161023, initial
%-

tc_track=[]; % init output

% poor man's version to check arguments
if ~exist('tc_track1','var'),tc_track1=[];end
if ~exist('tc_track2','var'),return;end 
if ~exist('start_year','var'),start_year=[];end 

% PARAMETERS
%

n_tracks1=length(tc_track1);
n_tracks2=length(tc_track2);

info.yyyy         =zeros(1,n_tracks1+n_tracks2);
info.datenum      =zeros(1,n_tracks1+n_tracks2);
info.track_i      =zeros(1,n_tracks1+n_tracks2);
info.track_source =zeros(1,n_tracks1+n_tracks2);

for track1_i=1:n_tracks1
    info.yyyy(track1_i)        =tc_track1(track1_i).yyyy(1);
    info.datenum(track1_i)     =tc_track1(track1_i).datenum(1);
    info.track_i(track1_i)     =track1_i;
    info.track_source(track1_i)=1;
end % track1_i

for track2_i=1:n_tracks2
    info.yyyy(n_tracks1+track2_i)        =tc_track2(track2_i).yyyy(1);
    info.datenum(n_tracks1+track2_i)     =tc_track2(track2_i).datenum(1);
    info.track_i(n_tracks1+track2_i)     =track2_i;
    info.track_source(n_tracks1+track2_i)=2;
end % track2_i

if isempty(start_year)
    start_year=min(info.yyyy);
elseif start_year==-1
    start_year=max(min(info.yyyy(1:n_tracks1)),min(info.yyyy(n_tracks1+1:end)));
end

% sort by datenum
[~,sort_index]   =sort(info.datenum);
info.datenum     =info.datenum(sort_index);
info.track_i     =info.track_i(sort_index);
info.track_source=info.track_source(sort_index);
info.yyyy        =info.yyyy(sort_index);


next_track=1;
tc_track=tc_track1(1);
end_year=-9999;
for track_i=1:length(info.datenum)
    
    if info.yyyy(track_i)>=start_year
        if info.track_source(track_i)==1
            tc_track(next_track)=tc_track1(info.track_i(track_i));
        else
            tc_track(next_track)=tc_track2(info.track_i(track_i));
        end
        end_year=max(end_year,tc_track(next_track).yyyy(1));
        next_track=next_track+1;
    end % year
end % track_i

fprintf('returning %i tracks since %i to %i\n',length(tc_track),start_year,end_year);

end % climada_tc_track_combine