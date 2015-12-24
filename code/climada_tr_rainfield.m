function [res, tc_track_ori, centroids] = climada_tr_rainfield(tc_track, centroids, equal_timestep, silent_mode, check_plot,unit_)
% TC rainfield calculation (rainsum)
% NAME:
%   climada_tc_rainfield
% PURPOSE:
%   given a single TC track (lat/lon,CentralPressure,MaxSustainedWind), calculate
%   the rain sum field at locations (=centroids)
%   mainly called from: see climada_tr_hazard_set
% CALLING SEQUENCE:
%   climada_tc_rainfield(tc_track, centroids)
% EXAMPLE:
%   res =
%   climada_tc_rainfield(tc_track(track_i),centroids)
%   climada_tc_rainfield(tc_track(track_i),centroids,1,1)
% INPUTS:
%   tc_track: a structure with the information for a single (1) tc track:
%       tc_track.lat
%       tc_track.lon
%       tc_track.MaxSustainedWind: maximum sustained wind speed (one-minute)
%       tc_track.MaxSustainedWindUnit as 'kn', 'mph', 'm/s' or 'km/h'
%       tc_track.CentralPressure: optional
%       tc_track.Celerity: translational (forward speed) of the hurricane.
%           optional, calculated from lat/lon if missing
%       tc_track.TimeStep: optional, only needed if Celerity needs to be
%           calculated, 6h assumed as default
%       tc_track.Azimuth: the forward moving angle, calculated if not given
%           to ensure consistency, it is even suggested not to pass Azimuth
%       tc_track.yyyy: 4-digit year, optional
%       tc_track.mm: month, optional
%       tc_track.dd: day, optional
%       tc_track.ID_no: unique ID, optional
%       tc_track.name: name, optional
%       tc_track.SaffSimp: Saffir-Simpson intensity, optional
%       NOTE: if empty, the user can also select a file with tc_track and will then
%       get promted for a single track number to use (mainly useful for TESTs)
%   centroids: a structure with the centroids information
%       centroids.lat: the latitude of the centroids
%       centroids.lon: the longitude of the centroids
% OPTIONAL INPUT PARAMETERS:
%   unit_: rainfall unit, either unit_='mm' (millimeters, default) or for inches unit_='in'
%   equal_timestep: if set=1 (default), first interpolate the track to a common
%       timestep, if set=0, no equalization of TC track data (not recommended)
%   silent_mode: if =1, do not write to stdout unless severe warning
%   check_plot: whether some plots are created =1 (default=0, silent_mode=1
%       sets check_plot also =0)
% OUTPUTS:
%   res.rainsum: the rain fall sum [mm per storm] at all centroids
%       the single-character variables refer to the Pioneer offering circular
%       that's why we kept these short names (so one can copy the OC for
%       documentation)
%   res.lat: the latitude of the centroids
%   res.lon: the longitude of the centroids
% RESTRICTIONS:
% MODIFICATION HISTORY:
% Lea Mueller, 20110606
% Martin Heyenn 20120503
% david.bresch@gmail.com, 20140804, GIT update
% david.bresch@gmail.com, 20141020, cleanup, inreach speedup
% David N. Bresch, david.bresch@gmail.com, 20150819, climada_global.centroids_dir introduced
%-

global climada_global
if ~climada_init_vars, return; end

if ~exist('tc_track'      ,'var'), tc_track       = []; end
if ~exist('centroids'     ,'var'), centroids      = []; end
if ~exist('equal_timestep','var'), equal_timestep = 1; end %don not set to []!
if ~exist('silent_mode'   ,'var'), silent_mode    = 1; end
if ~exist('check_plot'    ,'var'), check_plot     = 0; end
if ~exist('unit_'         ,'var'), unit_          = 'mm'; end

% PARAMETERS
%
% distance (in degree) around each node we process the rainfield
dlon=5; % default=5
dlat=dlon; % default=5

% prompt for tc_track if not given
if isempty(tc_track)
    tc_track = [climada_global.data_dir filesep 'tc_tracks' filesep '*.mat'];
    [filename, pathname] = uigetfile(tc_track, 'Select tc track set (usually a probabilistic one):');
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
    check_plot=1;
end



% prompt for centroids if not given
if isempty(centroids)
    centroids = [climada_global.centroids_dir filesep '*.mat'];
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

res          = []; % init output
tc_track_ori = tc_track;


if silent_mode, check_plot = 0; end

% refine tc_track to 1 hour timestep prior to windfield calculation
if equal_timestep
    if ~silent_mode,fprintf('NOTE: tc_track refined (1 hour timestep) prior to windfield calculation\n');end
    tc_track = climada_tc_equal_timestep(tc_track);
end

% calculate MaxSustainedWind if only CentralPressure given
if ~isfield(tc_track,'MaxSustainedWind') && isfield(tc_track,'CentralPressure')
    tc_track.MaxSustainedWind = tc_track.CentralPressure*0; % init
    tc_track.MaxSustainedWind(isnan(tc_track.MaxSustainedWind))=0;
end
% convert to kn
switch tc_track.MaxSustainedWindUnit
    case 'km/h'
        tc_track.MaxSustainedWind = tc_track.MaxSustainedWind/1.852;
    case 'mph'
        tc_track.MaxSustainedWind = tc_track.MaxSustainedWind/1.151;
    case 'm/s'
        tc_track.MaxSustainedWind = tc_track.MaxSustainedWind/1000*60*60/1.852;
    otherwise % already kn
end;
tc_track.MaxSustainedWindUnit='kn'; % after conversion



track_nodes_count = length(tc_track.lat);
centroid_count    = length(centroids.lat);
res.rainsum       = zeros(1,centroid_count);
res.lat           = centroids.lat';
res.lon           = centroids.lon';
if isfield(centroids,'OBJECTID')   ,res.OBJECTID = centroids.OBJECTID;   end
if isfield(centroids,'centroid_ID'),res.ID       = centroids.centroid_ID';end


% loop over track nodes
tic;
for i = 1:track_nodes_count
    
    % Define box where the CU's are in specific distance from the TC center
    %  inreach without on_land condition:
    
    % check for track nodes within centroids_rect (using inpolygon)
    %node_edges_x = [tc_track.lon(i)-dlon,tc_track.lon(i)-dlon,tc_track.lon(i)+dlon,tc_track.lon(i)+dlon,tc_track.lon(i)-dlon];
    %node_edges_y = [tc_track.lat(i)-dlat,tc_track.lat(i)+dlat,tc_track.lat(i)+dlat,tc_track.lat(i)-dlat,tc_track.lat(i)-dlat];
    %inreach = inpolygon(res.lon,res.lat,node_edges_x,node_edges_y);
    
    % check for track nodes within centroids_rect (using many logical)
    %inreach = res.lon > (tc_track.lon(i) - 5) & ...
    %      res.lon < (tc_track.lon(i) + 5) & ...
    %      res.lat > (tc_track.lat(i) - 5) & ...
    %      res.lat < (tc_track.lat(i) + 5);
    
    inreach=abs(res.lon-tc_track.lon(i))<dlon & abs(res.lat-tc_track.lat(i))<dlat;
    
    % Calculate distance for individual gridpoints from TC Center for windfield
    %  calculation/orographic rain & R-CLIPER - only calculated for the
    %  inreach-box, outside set to zero.
    if any(inreach)
        [fRadius_km GridVect] = climada_nonspheric_distance_m(...
            res.lon(inreach),...
            res.lat(inreach),...
            tc_track.lon(i),...
            tc_track.lat(i),...
            centroid_count,...
            inreach);
        
        % calculate rain rate in mm / h
        rainrate         = climada_RCLIPER(tc_track.MaxSustainedWind(i),...
            inreach,...
            fRadius_km);
        res.rainrate(i,:)= sparse(rainrate); % uses a LOT of time, unnecessary
        res.rainsum      = res.rainsum + rainrate;  %total sum of mm per wind storm
    end
    
end %Loop over all TC Nodes


title_str = [tc_track.name ', ' datestr(tc_track.datenum(1))];
if ~silent_mode, fprintf('%f secs for %s rainfall sum field\n',toc,deblank(title_str));end


%--------------
% FIGURE
%--------------
if check_plot
    fprintf('preparing rainfall sum footprint plot\n')
    
    
    %scale figure according to range of longitude and latitude
    scale  = max(centroids.lon) - min(centroids.lon);
    scale2 =(max(centroids.lon) - min(centroids.lon))/...
        (min(max(centroids.lat),60)-max(min(centroids.lat),-50));
    height = 0.5;
    if height*scale2 > 1.2; height = 1.2/scale2; end
    fig = climada_figuresize(height,height*scale2+0.15);
    
    
    %check for unit input
    if strcmp(unit_,'mm')
        unit=1  ;            %mm
        unitstr=' mm';
    elseif strcmp(unit_,'in')
        unit=0.0393700787; %inch
        unitstr=' in';
    end
    
    
    % create gridded values
    [X, Y, gridded_VALUE] = climada_gridded_VALUE(res.rainsum, centroids);
    gridded_max           = max(max(gridded_VALUE)*unit);
    %MH set vlaues lower than 0.1*unit to NaN for ploting
    gridded_VALUE(gridded_VALUE<(0.1*unit)) = NaN;
    
    %plot
    contourf(X, Y, full(gridded_VALUE)*unit,50,'edgecolor','none') %mm to in *0.0393700787
    hold on
    climada_plot_world_borders(0.7)
    climada_plot_tc_track_stormcategory(tc_track_ori);
    
    
    %plot centroids?
    plot(centroids.lon, centroids.lat, '+r','MarkerSize',0.8,'linewidth',0.1)
    
    axis equal
    axis([min(centroids.lon)-scale/30  max(centroids.lon)+scale/30 ...
        max(min(centroids.lat),-50)-scale/30  min(max(centroids.lat),60)+scale/30])
    
    
    if log10(gridded_max)>2
        caxismax=(floor(gridded_max/100)+1)*100;
        caxis([0 caxismax])
    else
        caxismax=(floor(gridded_max/10)+1)*10;
        caxis([0 caxismax])
    end
    
    %caxis([0 600])
    
    
    %colorbar
    
    %number of colors
    %steps10=(floor(gridded_max/100)+1)*100/10;
    steps10=10;
    %steps10=8;
    
    
    %create colormap
    cmap1=[];
    cmap2=[];
    
    startcolor  =[0.89	0.93	0.89];
    middlecolor1=[0.55	0.78	0.59];
    middlecolor2=[0.43	0.84	0.78];
    endcolor    =[0.05	0.37	0.55];
    
    for i=1:3
        cmap1(:,i)=startcolor(i):(middlecolor1(i)-startcolor(i))/(steps10/2-1):middlecolor1(i);
        cmap2(:,i)=middlecolor2(i):(endcolor(i)-middlecolor2(i))/(steps10/2-1):endcolor(i);
    end
    
    cmap=[cmap1;cmap2];
    colormap(cmap)
    
    c=colorbar('YTick',0:caxismax/10:caxismax);
    %c=colorbar('YTick',0:75:600);
    ylabel(c,['rainsum [',unitstr,']'],'fontsize',8)
    set(gca,'fontsize',8)
    
    %labels
    xlabel('Longitude','fontsize',8)
    ylabel('Latitude','fontsize',8)
    
    %title(title_str,'interpreter','none')
    stormdate     = tc_track.datenum(1);
    stormduration = tc_track.datenum(end) - tc_track.datenum(1);
    stormname     = tc_track.name;
    stormname(stormname == '_') = ' ';
    
    %find centroid with maximum rain sum and derive duration of rain
    [rainsum_max centroid_max] = max(res.rainsum);
    pos                        = find(res.rainrate(:,centroid_max)>0);
    pos_count                  = length(pos); %No. of hours with rain
    title({[stormname ', ' datestr(stormdate,'dd mmm yyyy') ', duration ' datestr(stormduration,'dd HH') ' (days hours)'];...
        ['max:       ' int2str(round(gridded_max)), unitstr,' in ' datestr(pos_count/24,'dd HH') ' (days hours)'];...
        ['average:   ' num2str(gridded_max/pos_count,'%10.1f'),unitstr,' h^{-1}                           ']},'fontsize',8)
    
    choice = questdlg('print?','print');
    switch choice
        case 'Yes'
            check_printplot = 1;
        case 'No'
            check_printplot = 0;
        case 'Cancel'
            return
    end
    
    if check_printplot
        %foldername = ['\results\mozambique\V\footprint_' int2str(storm_no) '_' storm_name '.pdf'];
        foldername = [filesep 'results' filesep 'footprint_rainsum_' tc_track.name '.pdf'];
        print(fig,'-dpdf',[climada_global.data_dir foldername])
        %close
        fprintf('saved 1 FIGURE in folder %s \n', foldername);
    end
end


return