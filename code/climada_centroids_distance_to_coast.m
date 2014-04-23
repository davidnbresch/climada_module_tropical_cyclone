function centroids = climada_centroids_distance_to_coast(centroids, coastline, check_figure)

% calculate distance to coast in km for every centroid and add information
% in the field centroids.dist_to_coast
% NAME:
%   climada_centroids_distance_to_coast
% PURPOSE:
%   calculate distance to coast in km for every centroid and add information
%   in the field centroids.dist_to_coast
%   within: climada_hazard_distance_to_coast
% CALLING SEQUENCE:
%   centroids = climada_centroids_distance_to_coast(centroids, coastline,
%   check_figure)
% EXAMPLE:
%   centroids = climada_centroids_distance_to_coast
% INPUTS:
%   none, if coastline_file empty default file is loaded from 
% OPTIONAL INPUT PARAMETERS:
%   check_figure to create figure
% OUTPUTS:
%   centroids structure variable including field centroids.dist_to_coast
%   with distance to coast in km
% RESTRICTIONS:
% MODIFICATION HISTORY:
% Lea Mueller, 20121203
%-

% init global variables
global climada_global
if ~climada_init_vars,return;end

% check inputs
if ~exist('centroids'   , 'var'), centroids    = []; end
if ~exist('coastline'   , 'var'), coastline    = []; end
if ~exist('check_figure', 'var'), check_figure = []; end


if isempty(coastline)
    coastline = climada_coastline_read;
end

if isempty(centroids)
    centroids            = [climada_global.system_dir filesep '*.mat'];
    centroids_default    = [climada_global.system_dir filesep 'Select centroids .mat'];
    [filename, pathname] = uigetfile(centroids, 'Select tc track set:',centroids_default);
    if isequal(filename,0) || isequal(pathname,0)
        return; % cancel
    else
        centroids = fullfile(pathname,filename);
    end
end
if ~isstruct(centroids) % load, if filename given
    centroids_file = centroids; centroids = [];
    vars = whos('-file', centroids_file);
    load(centroids_file);
    if ~strcmp(vars.name,'centroids')
        centroids = eval(vars.name);
        clear (vars.name)
    end
end

[a, b] = size(centroids.Longitude);

if isfield(centroids,'onLand')
    onLandindex = centroids.onLand == 1;
else
    onLandindex = ones([a b]);
end
index_ = find(onLandindex);
    
if isfield(centroids,'dist_to_coast')
    fprintf('Distance to coast is already calculated for these centroids.\n')
else
    % initialize distance to coast with zeros
    
    centroids.dist_to_coast = zeros([a b]);
    min_dd                  = zeros([a b]);

    no_centroids = length(centroids.Longitude(onLandindex));
    t0           = clock;
    msgstr       = sprintf('Add distance to coast for %d centroids\n', no_centroids);
    h            = waitbar(0, msgstr, 'Name','Calculate distance to coast for centroids');
    mod_step     = 10; % first time estimate after 10 tracks, then every 100

    for cen_i = 1:no_centroids
        %calculate distance to all coastline points and find closest, and
        %therrefore minimal distance to coast (in km)
        cen_ii = index_(cen_i);
        dd     = climada_geo_distance(centroids.Longitude(cen_ii), centroids.Latitude(cen_ii),...
                                      coastline.lon, coastline.lat)/1000;                          
        [min_dd(cen_i), pos] = min(dd);

        if mod(cen_i,mod_step)==0
            mod_step              = 100;
            t_elapsed_centroids   = etime(clock,t0)/cen_i;
            centroids_remaining   = no_centroids - cen_i;
            t_projected_centroids = t_elapsed_centroids * centroids_remaining;
            msgstr                = sprintf('est. %i seconds left (%i/%i centroids)',ceil(t_projected_centroids), cen_i, no_centroids);
            waitbar(cen_i/no_centroids,h,msgstr); % update waitbar
        end
    end
    centroids.dist_to_coast(logical(onLandindex)) = min_dd;
    close(h); % dispose waitbar
end

if check_figure
    delta   = 2;
    axislim = [min(centroids.Longitude(onLandindex))-delta  max(centroids.Longitude(onLandindex))+delta ...
               min(centroids.Latitude (onLandindex))-delta  max(centroids.Latitude (onLandindex))+delta];
    fig_relation = diff(axislim(1:2))/diff(axislim(3:4));
    if fig_relation<1
        fig_height = 0.6; 
        fig_width  = fig_height*fig_relation; 
    else
        fig_width  = 0.8;
        fig_height = fig_width/fig_relation; 
    end   
    fig = climada_figuresize(fig_height, fig_width);
    climada_plot_world_borders
    cmap = flipud(jet);
    cbar = plotclr(centroids.Longitude(logical(onLandindex)), centroids.Latitude(logical(onLandindex)), centroids.dist_to_coast(logical(onLandindex)),...
                   's', 2, 1, [], [], cmap, [], 1);
    set(get(cbar,'ylabel'),'String', 'Distance to coast (km)','fontsize',11)
    hold on
    plot(centroids.Longitude(~onLandindex), centroids.Latitude(~onLandindex),'.','color',[139 136 120]/255)
    axis equal
    axis(axislim)
end

                          
                          
                          
                          
                          
                          
                          

