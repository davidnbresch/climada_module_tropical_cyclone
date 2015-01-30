
function climada_tc_track_windfield_google_earth(tc_track,centroids,aggregation,google_earth_save)
% climada_tc_track_windfield_google_earth
% MODULE:
%   tc_hazard_advanced
% NAME:
%   climada_tc_track_windfield_google_earth
% PURPOSE:
%   create kmz-file of wind field animation in google earth
% CALLING SEQUENCE:
%   climada_tc_track_windfield_google_earth(tc_track, centroids, aggregation, google_earth_save)
% EXAMPLE:
%   climada_tc_track_windfield_google_earth
% INPUTS:
%   tc_track: track set, historical or probabilistic, but only historical
%   tracks are visualized. prompted for, if not given.
%   centroids
% OPTIONAL INPUT PARAMETERS:
%   aggregation:        aggregation time, default is 6h
%   google_earth_save:  filename, to be saved in
%                       \climada\data\tc_tracks\tc_track....kml
% OUTPUTS:
%   kmz-file, visualisation of historical track and winfield with time stamp in google
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


%% prompt for google_earth_save file
if isempty(google_earth_save) % local GUI
    google_earth_save = [climada_global.data_dir filesep 'tc_tracks' filesep 'Select name to save google earth visualiation tc_tracks_wind field animation.kmz'];
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
k = kml(google_earth_save);


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


%% colormap, constant color range
[cmap c_ax]= climada_colormap('TC');
cmap_ori = cmap;
c_ax_ori = c_ax;
% do show wind speeds/intensities only after reaching a certain level
if c_ax(1)==0
    levels = linspace(c_ax(1), c_ax(2), length(cmap)+1);
    c_ax(1) = levels(2);
    cmap(1,:)=[];
end

for agg_i = 1:aggregation_count
    [X, Y, gridded_VALUE] = climada_gridded_VALUE(full(res.gust_aggr(agg_i,:)), centroids);
    if sum(gridded_VALUE(:)>10)>0    
        
        % to debug
        %figure
        %contourf(X, Y, gridded_VALUE, levels(2:end))
        %colormap(cmap)
        %colorbar
        %title(datestr(tc_track.datenum(1)+agg_i*aggregation/24,'dd mmm yyyy, HHpm'))
        %caxis(c_ax)
        %hold on
        %climada_plot_world_borders([],[],[],1)
        
        kk = k.createFolder(datestr(tc_track.datenum(1)+agg_i*aggregation/24,'dd mmm yyyy, HHpm'));
        
        kk.contourf(X,Y,gridded_VALUE,...
             'name',tc_track.name,'description','test',...
             'colorMap',cmap,'lineColor','00FFFFFF','transparency',0.5,...
             'caxis',c_ax,...
             'timeSpanBegin',datestr(stormdate+(agg_i-1)*aggregation/24, 'yyyy-mm-ddTHH:MM:SSZ'),...
             'timeSpanEnd',datestr(stormdate+(agg_i  )*aggregation/24, 'yyyy-mm-ddTHH:MM:SSZ'));
    end
end


%% visualize track lines and nodes
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
kk = k.createFolder(tc_track.name);
description_str = sprintf('%s, %s - %s', tc_track.name, ...
                        datestr(tc_track.datenum(1),'dd mmm yyyy'), ...
                        datestr(tc_track.datenum(end),'dd mmm yyyy'));
for node_i = 1:length(tc_track.lon)-1
    v       = tc_track.MaxSustainedWind(node_i);
    v_color = find (v < v_categories);
    v_color = v_color(1);


    kk.plot(tc_track.lon(node_i:node_i+1), tc_track.lat(node_i:node_i+1),...
         'name',tc_track.name,'description',description_str,...
         'lineColor',colors_(v_color,:),...
         'timeSpanBegin',datestr(tc_track.datenum(node_i),'yyyy-mm-ddTHH:MM:SSZ'),...
         'timeSpanEnd' ,datestr(tc_track.datenum(end)+6/24,'yyyy-mm-ddTHH:MM:SSZ'));

    kk.point(tc_track.lon(node_i), tc_track.lat(node_i),100,...
         'description',description_str,...
         'iconURL','http://maps.google.com/mapfiles/kml/shapes/donut.png',...
         'iconScale',0.5,...
         'iconColor',colors_(v_color,:),...
         'timeSpanBegin',datestr(tc_track.datenum(node_i),'yyyy-mm-ddTHH:MM:SSZ'),...
         'timeSpanEnd' ,datestr(tc_track.datenum(end)+6/24,'yyyy-mm-ddTHH:MM:SSZ')); 
end

kk.point(tc_track.lon(node_i+1), tc_track.lat(node_i+1),100,...
        'description',description_str,...
        'iconURL','http://maps.google.com/mapfiles/kml/shapes/donut.png',...
        'iconScale',0.5,...
        'iconColor',colors_(v_color,:),...
        'timeSpanBegin',datestr(tc_track.datenum(node_i),'yyyy-mm-ddTHH:MM:SSZ'),...
        'timeSpanEnd' ,datestr(tc_track.datenum(end)+6/24,'yyyy-mm-ddTHH:MM:SSZ')); 

    
%% open visualition in google earth
k.run
   





