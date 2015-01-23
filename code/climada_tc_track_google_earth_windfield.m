
function climada_tc_track_google_earth_windfield(tc_track,centroids,aggregation,google_earth_save)

% ----------------------------------------
%% create kml-file of wind field animation
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
%   google_earth_save:  filename, to be saved in
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
if ~exist('centroids'        ,'var'), centroids         = []; end
if ~exist('aggregation'      ,'var'), aggregation       = []; end
if ~exist('google_earth_save','var'), google_earth_save = []; end
if isempty(aggregation)             , aggregation       = 6 ; end


%% prompt for tc_track if not given
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
    prompt   ='Type specific No. of track to print windfield [e.g. 1, 10, 34]:';
    name     =' No. of track';
    defaultanswer = {'34'};
    answer = inputdlg(prompt,name,1,defaultanswer);
    track_no = str2num(answer{1});
    tc_track = tc_track(track_no);
end

%% prompt for centroids if not given
if isempty(centroids)
    centroids = [climada_global.system_dir filesep '*.mat'];
    [filename, pathname] = uigetfile(centroids, 'Select centroids:');
    if isequal(filename,0) || isequal(pathname,0)
        return; % cancel
    else
        centroids = fullfile(pathname,filename);
    end
end
% load the centroids, if a filename has been passed
if ~isstruct(centroids)
    centroids_file = centroids;
    centroids      = [];
    vars = whos('-file', centroids_file);
    load(centroids_file);
    if ~strcmp(vars.name,'centroids')
        centroids = eval(vars.name);
        clear (vars.name)
    end
end

                                                             
%%

% colormap, constant color range
load ([climada_global.system_dir filesep 'colormap_gray_red_90'])
colormap(gray_red)
max_wind = 90;
caxis([0 max_wind])

% colors according to saffir-simpson scale
v_categories = [34 64 83 96 113 135 1000];
colors_      = [ 'ffffaa00' ; %blue
                 'ff00aa55' ; %green
                 'ff00ffff' ; %yellow
                 'ff00aaff' ; %dark yellow  
                 'ff0055ff' ; %orange  
                 'ff0000ff' ; %red
                 'ff000079' ; %dark red
                ];
            
kmlStr  = [];
kmlStr_ = [];
kmlStr__= [];
values = [1:10:max_wind];


%---------------------------
%% Calculations
%---------------------------
% wind field for every hour (for every node from tc_track)
% equal timestep within this routine
res = climada_tc_windfield_timestep(tc_track,centroids,1); 
stormdate = tc_track.datenum(1);
stormname = tc_track.name;
stormname(stormname == '_') = ' ';

% aggregate wind field for specific hours (unit remains mm/s)
[a b]             = size(res.gust);
aggregation_count = floor(a/aggregation);
if aggregation > 1
    for i = 1:aggregation_count
        res.gust_aggr(i,:) = mean(res.gust((i-1)*aggregation+1:i*aggregation,:));
    end
else
    res.gust_aggr = res.gust;
end


for agg_i = 1:aggregation_count
    [X, Y, gridded_VALUE] = climada_gridded_VALUE(full(res.gust_aggr(agg_i,:)), centroids);
    if sum(gridded_VALUE(:)>10)>0         
        kmlStr = [kmlStr ...
                  ge_contourf(X,Y,gridded_VALUE,...
                        'cMap',gray_red,...
                        'lineValues',values,...
                        'lineColor','00FFFFFF',...
                        'timeSpanStart',datestr(stormdate+(agg_i-1)*aggregation/24, 'yyyy-mm-ddTHH:MM:SSZ'),...
                        'timeSpanStop' ,datestr(stormdate+(agg_i  )*aggregation/24, 'yyyy-mm-ddTHH:MM:SSZ'),...
                        'name',stormname,...
                        'description',[stormname int2str(i) ' highest'],...
                        'polyAlpha','6f'...
                        )];

    end
end

for node_i = 1:length(tc_track.lon)-1
        v       = tc_track.MaxSustainedWind(node_i);
        v_color = find (v < v_categories);
        v_color = v_color(1);
        kmlStr_ = [kmlStr_ ...
                  ge_plot(tc_track.lon(node_i:node_i+1), tc_track.lat(node_i:node_i+1),...
                         'name',tc_track.name,...
                         'lineColor',colors_(v_color,:),...
                         'timeSpanStart',datestr(tc_track.datenum(node_i),'yyyy-mm-ddTHH:MM:SSZ'),...
                         'timeSpanStop' ,datestr(tc_track.datenum(end)+6/24,'yyyy-mm-ddTHH:MM:SSZ')...
                         )];
end
    
for node_i = 1:length(tc_track.lon)
        v       = tc_track.MaxSustainedWind(node_i);
        v_color = find (v < v_categories);
        v_color = v_color(1);
        kmlStr__= [kmlStr__...
              ge_point(tc_track.lon(node_i), tc_track.lat(node_i),100,...
              'iconURL','http://maps.google.com/mapfiles/kml/shapes/donut.png',...
              'iconColor',colors_(v_color,:),...
              'timeSpanStart',datestr(tc_track.datenum(node_i),'yyyy-mm-ddTHH:MM:SSZ'),...
              'timeSpanStop' ,datestr(tc_track.datenum(end)+6/24,'yyyy-mm-ddTHH:MM:SSZ'))];
end

    
    
%% prompt for google_earth_save file
if isempty(google_earth_save) % local GUI
    google_earth_save = [climada_global.data_dir filesep 'tc_tracks' filesep 'Select name to save google earth visualiation tc_tracks_wind field animation.kml'];
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

ge_output(google_earth_save,[kmlStr kmlStr_ kmlStr__])





%%



