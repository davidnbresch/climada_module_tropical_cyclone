function climada_hazard_footprint_plot(hazard, event_no, nametag)
% TC footprint figure
% NAME:
%   climada_plot_tc_footprint
% PURPOSE:
% create footprint figure
% CALLING SEQUENCE:
%   [contr t_handle] = climada_plot_tc_footprint(hazard, tc_track, track_no)
% EXAMPLE:
%   climada_plot_tc_footprint
% INPUTS:
%   hazard: hazard.intensity with intensities per centroid
%   event_no: a number, event number of hazard to be shown
% OPTIONAL INPUT PARAMETERS:
% OUTPUTS:
%   figure with footprint
% RESTRICTIONS:
% MODIFICATION HISTORY:
% Lea Mueller, muellele@gmail.com, 20121205
% Lea Mueller, muellele@gmail.com, 20150423
%-

global climada_global
if ~climada_init_vars, return; end
if ~exist('hazard'       , 'var'), hazard       = []; end
if ~exist('event_no'     , 'var'), event_no     = []; end
if ~exist('nametag'      , 'var'), nametag      = ''; end

% prompt for hazard if not given
if isempty(hazard)
    hazard               = [climada_global.data_dir filesep 'hazard' filesep '*.mat'];
    hazard_default       = [climada_global.data_dir filesep 'hazard' filesep 'Select hazard .mat'];
    [filename, pathname] = uigetfile(hazard, 'Select hazard set:',hazard_default);
    if isequal(filename,0) || isequal(pathname,0)
        return; % cancel
    else
        hazard = fullfile(pathname,filename);
    end
end
% load the hazard set, if a filename has been passed
if ~isstruct(hazard)
    hazard_file = hazard;
    hazard      = [];
    vars = whos('-file', hazard_file);
    load(hazard_file);
    hazard=climada_hazard2octave(hazard); % Octave compatibility for -v7.3 mat-files
    if ~strcmp(vars.name,'hazard')
        hazard = eval(vars.name);
        clear (vars.name)
    end
    prompt   ='Type specific No. of track to print windfield [e.g. 1, 10, 34, 1011]:';
    name     =' No. of track';
    defaultanswer = {'1'};
    answer   = inputdlg(prompt,name,1,defaultanswer);
    event_no = str2double(answer{1});
    intensity_field = hazard.intensity(event_no,:);
end

if isempty(event_no); event_no = 1; end
if event_no>size(hazard.intensity,1); event_no = size(hazard.intensity,1);end
intensity_field = hazard.intensity(1,:);

% some hazard information that will be checked
if ~isfield(hazard,'peril_ID'); hazard.peril_ID = ''; end
if ~isfield(hazard,'units'    ); hazard.units   = ''; end

% set the axis
% longitude, latitude range
lon     = hazard.lon;
lat     = hazard.lat;
dlon    = max(lon) - min(lon);
dlat    = max(lat) - min(lat);
x_range = [min(lon)-dlon/10 max(lon)+dlon/10];
y_range = [min(lat)-dlat/10 max(lat)+dlat/10];
hb      = get(gca,'PlotBoxAspectRatio');
hb      = hb(1)/hb(2);
if hb/(diff(x_range)/diff(y_range))<1
    dif     = ( diff(x_range)/hb-diff(y_range) )/2;
    y_range = [y_range(1)-dif y_range(2)+dif];
    set(gca,'xlim',x_range,'ylim',y_range)
else
    dif     = ( diff(y_range)*hb-diff(x_range) )/2;
    x_range = [x_range(1)-dif x_range(2)+dif]; 
    set(gca,'xlim',x_range,'ylim',y_range)
end

% plot figure only if there exist intensity values
intesity_cutoff   = 0.001;
gridded_max_round = [];
%special settings for tc wind 
if strcmp(hazard.peril_ID,'TC'); 
    intesity_cutoff = 5; 
    gridded_max_round = 70; %45
end 
if size(find(intensity_field>intesity_cutoff),2)>intesity_cutoff %any(full(windfield))  
    % create colormap according to peril ID
    [cmap c_ax] = climada_colormap(hazard.peril_ID);
    % cmap_= [1.0000    1.0000    1.0000;        
    %         0.8100    0.8100    0.8100; %1.0000    1.0000    1.0000;            
    %         0.6300    0.6300    0.6300;
    %         1.0000    0.8000    0.2000;
    %         0.9420    0.6667    0.1600; %0.8839    0.5333    0.1200;
    %         0.8259    0.4000    0.0800;
    %         0.7678    0.2667    0.0400;
    %         0.7098    0.1333         0;
    %         0.5020         0    0.5020;
    %         0.2941         0    0.5098];
    %cbar = plotclr(centroids.lon, centroids.lat, full(windfield), 's',4,1,-1,[],cmap_,0,0);
    %set(get(cbar,'ylabel'),'String', 'Wind speed (m s^{-1})' ,'fontsize',8);
    
    % % create gridded values
    intensity_field = full(intensity_field);
    w_index         = intensity_field>intesity_cutoff;
    centroids.lon   = hazard.lon(w_index);
    centroids.lat   = hazard.lat(w_index);
    [X, Y, gridded_VALUE] = climada_gridded_VALUE(intensity_field(w_index), centroids);
    gridded_VALUE(isnan(gridded_VALUE)) = 0;
    gridded_max       = max(max(gridded_VALUE));
    if isempty(gridded_max_round); gridded_max_round = gridded_max; end %45; 
    %bins = gridded_max_round/size(cmap,1)*1.2;
    bins = gridded_max_round/size(cmap,1);
    hold on
    contourf(X, Y, full(gridded_VALUE),...
                    0:bins:gridded_max_round,'edgecolor','none')
    caxis([0 gridded_max_round]) 
    colormap(cmap)
    %plot(centroids.lon, centroids.lat,'+')
    %colorbartick           = [0:10:gridded_max_round round(gridded_max)];
    % colorbartick           = [0:bins:gridded_max_round round(gridded_max)];
    % colorbarticklabel      = num2cell(colorbartick);
    % colorbarticklabel{end} = [num2str(gridded_max,'%10.2f') 'max'];
    % colorbarticklabel{end} = [int2str(gridded_max)          'max'];
    % t = colorbar('YTick',colorbartick,'yticklabel',colorbarticklabel);
    t = colorbar;
    cbar_label = sprintf('Intensity %s (%s)', hazard.peril_ID, hazard.units);
    set(get(t,'ylabel'),'String', cbar_label,'fontsize',10);
end
climada_plot_world_borders(0.7,[], [], 1)


% climada_plot_tc_track_stormcategory(tc_track, 4, [], 1);
% % plot(centroids.lon,centroids.lat,'.k','markersize',6,'linewidth',1)
% 
% if ~isfield(tc_track,'category'),tc_track=climada_tc_stormcategory(tc_track);end
% 
% title_str = sprintf('tc track %s: %s, \t %s - %s, category %d', nametag, tc_track.name, ...
%                     datestr(tc_track.datenum(1),'dd/mmm') ,datestr(tc_track.datenum(end),'dd/mmm/yyyy'),...
%                     tc_track.category);               
% title(title_str,'interpreter','none','fontsize',8)


% if isempty(check_printplot)
%     choice = questdlg('print?','print');
%     switch choice
%     case 'Yes'
%         check_printplot = 1;
%     case 'No'
%         check_printplot = 0;
%     case 'Cancel'
%         return
%     end
% end
% 
% if check_printplot %(>=1)   
%     foldername = [filesep 'results' filesep 'footprint_' tc_track.name '.pdf'];
%     print(fig,'-dpdf',[climada_global.data_dir foldername])
%     %close
%     fprintf('saved 1 FIGURE in folder %s \n', foldername);
% end

 
end
    
