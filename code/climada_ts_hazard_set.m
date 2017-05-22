function hazard=climada_ts_hazard_set(hazard,hazard_set_file,elevation_data,check_plots)
% climada storm surge TS hazard event set
% NAME:
%   climada_ts_hazard_set
% PURPOSE:
%   create a storm surge (TS) hazard event, based on an existing
%   tropical cyclone (TC) hazard event set
%
%   two steps:
%   1) check wether we need to obtain bathymetry (high res, calls etopo_get
%   from module etopo) or elevation is provided (see elevation_data)
%   2) convert all TC footprints into TS footprints
%   see CORE_CONVERSION in code below for the conversion formula
%
%   Note on bathymetry data: see https://github.com/davidnbresch/climada_module_elevation_models 
%
%   see tc_surge_TEST for a testbed for this code
%
% CALLING SEQUENCE:
%   hazard=climada_ts_hazard_set(hazard,hazard_set_file,elevation_data,check_plot)
% EXAMPLE:
%   hazard_TC=climada_hazard_load('TCNA_today_small')
%   hazard_TS=climada_ts_hazard_set(hazard_TC,'TCNA_today_small_TS')
% INPUTS:
%   hazard: an already existing tropical cyclone (TC) hazard event set (a
%       TC hazard structure)
%       > prompted for if not given (for .mat file containing a TC hazard event set)
%       Note: if hazard.elevation_m exists on input, this elevation
%       information is used, hence elevation_data is ignored
%       The variable hazard is modified on output (saves a lot of memory).
%   hazard_set_file: the name of the newly created storm surge (TS) hazard
%       event set (if ='NO_SAVE', the hazard is just returned, not saved)
%       > promted for if not given
% OPTIONAL INPUT PARAMETERS:
%   elevation_data: if a scalar, take elevation (or bathymetry) from etopo
%       (needs module etopo). If =1, save bathymetry in a .mat file (speeds up
%       subsequent calls), if =0 do not save bathymetry (default)
%       If elevation_data is a structure with the fields lon(i), lat(i),
%       elevation_m(i) and centroid_ID(i), check for centroid IDs being the
%       same as in hazard.centroid_ID and then just 'attach' the elevation
%       to hazard (i.e. hazard.elevation_m=elevation_data.elevation_m)
%       Note: if hazard.elevation_m exists on input, this elevation
%       information is used.
%   check_plots: =1, do show check plots (only if BATI used), =0: no plots (default)
% OUTPUTS:
%   hazard: a hazard event set, see core climada doc
%       also written to a .mat file (see hazard_set_file)
%       NOTE: for memory allocation reasons, the input hazard is used and
%       modified to create the output hazard
% MODIFICATION HISTORY:
% david.bresch@gmail.com, 20140421
% david.bresch@gmail.com, 20141017, module path relative
% david.bresch@gmail.com, 20141026, save_bathymetry_flag
% david.bresch@gmail.com, 20150106, elevation_data
% david.bresch@gmail.com, 20160516, elevation_data single precision allowed (as eg from SRTM)
% david.bresch@gmail.com, 20160525, better error messaging for ETOPO issues
% david.bresch@gmail.com, 20160529, renamed to climada_ts_hazard_set and tc_surge_hazard_create deleted
% david.bresch@gmail.com, 20161009, strcmpi used
%-

global climada_global
if ~climada_init_vars,return;end % init/import global variables

% poor man's version to check arguments
if ~exist('hazard','var'),hazard=[];end
if ~exist('hazard_set_file','var'),hazard_set_file='';end
if ~exist('elevation_data','var'),elevation_data=[];end
if ~exist('check_plots','var'),check_plots=0;end

% PARAMETERS
%
module_data_dir=[fileparts(fileparts(mfilename('fullpath'))) filesep 'data'];
if ~isdir(module_data_dir),mkdir(fileparts(module_data_dir),'data');end % create the data dir, should it not exist (no further checking)
%
% whether we save the bathymetry tile for subsequent use
save_bathymetry_flag=0; % default=0, see elevation_data

% prompt for TC hazard event set if not given
if isempty(hazard) % local GUI
    TC_hazard_set_file=[climada_global.data_dir filesep 'hazards' filesep '*.mat'];
    [filename, pathname] = uigetfile(TC_hazard_set_file, 'Select a TC hazard event set:');
    if isequal(filename,0) || isequal(pathname,0)
        return; % cancel
    else
        TC_hazard_set_file=fullfile(pathname,filename);
    end
    load(TC_hazard_set_file);
end

% prompt for hazard_set_file if not given
if isempty(hazard_set_file) % local GUI
    hazard_set_file=[climada_global.data_dir filesep 'hazards' filesep 'TS_hazard.mat'];
    [filename, pathname] = uiputfile(hazard_set_file, 'Save new TS hazard event set as:');
    if isequal(filename,0) || isequal(pathname,0)
        return; % cancel
    else
        hazard_set_file=fullfile(pathname,filename);
    end
end

% complete path, if missing
[fP,fN,fE]=fileparts(hazard_set_file);
if isempty(fP),hazard_set_file=[climada_global.data_dir filesep 'hazards' filesep fN fE];end

if ~isstruct(elevation_data)
    save_bathymetry_flag=elevation_data;
    elevation_data=[]; % make empty
end

% 1) try to obtain the elevation in different ways
% ------------------------------------------------

if ~isfield(hazard,'elevation_m') && ~isempty(elevation_data)
    % first, use data as provided in elevation_data
    if isfield(elevation_data,'centroid_ID')
        if length(elevation_data.centroid_ID)==length(hazard.centroid_ID)
            if sum(elevation_data.centroid_ID-hazard.centroid_ID)==0
                % elevation provided at hazard centroids
                hazard.elevation_m=elevation_data.elevation_m;
            end
        end
    end
end

if ~isfield(hazard,'elevation_m')
    
    if ~isempty(elevation_data)
        % if elevation_data provided, we should not get here
        fprintf('Warning: elevation_data failed, using etopo\n');
    end
    
    % second, if elevation_data is empty, use the elevation module

    if isempty(which('etopo_get')) % check for elevation module
        fprintf(['WARNING: install climada elevation module first. Please download ' ...
            '<a href="https://github.com/davidnbresch/climada_module_elevation_models">'...
            'climada_module_elevation_models</a> from Github.\n'])
        fprintf('Note: code %s continues, but might encounter problems\n',mfilename);
    end
        
    % prep the region we need (rectangular region encompassing the hazard centroids)
    centroids_rect=[min(hazard.lon) max(hazard.lon) min(hazard.lat) max(hazard.lat)];
    
    % 1) create the bathymetry file
    % -----------------------------
    
    bb=1; % 1 degree of bathy outside centroids (to allow for smooth interp)
    bathy_coords=[centroids_rect(1)-bb centroids_rect(2)+bb centroids_rect(3)-bb centroids_rect(4)+bb];
    
    % cut the bathymtery data out of the global topography dataset
    
    % file the bathymetry gets stored in ( might not be used later on, but
    % saved to speed up re-creation of the TS hazard event set
    [~,fN]=fileparts(hazard_set_file);
    Bathymetry_file=[module_data_dir filesep fN '_bathy.mat'];
    
    if ~exist(Bathymetry_file,'file')
        
        if ~exist('etopo_get','file')
            % safety to inform the user in case he misses the ETOPO module
            fprintf(['ERROR: etopo_get function not found, install climada elevation module first. Please download ' ...
                '<a href="https://github.com/davidnbresch/climada_module_elevation_models">'...
                'climada_module_elevation_models</a> from Github.\n'])
            hazard=[];
            return
        end
        BATI=etopo_get(bathy_coords);
        if isempty(BATI),hazard=[];
            return
        end % error messages from etopo_get already
        
        if save_bathymetry_flag
            fprintf('saving bathymetry as %s (you might later delete this file)\n',Bathymetry_file);
            save(Bathymetry_file,'BATI');
        end
        
        if check_plots
            figure('Name','Bathymetry','Color',[1 1 1]);
            pcolor(BATI.x,BATI.y,BATI.h)
            hold on
            shading interp
            colorbar
            caxis([-200 20]) % a reasonable colorbar for elevation
            axis equal
            climada_plot_world_borders
            axis(bathy_coords);
        end
        
    else
        fprintf('reading bathymetry from %s\n',Bathymetry_file);
        load(Bathymetry_file); % contains BATI
    end
    
    hazard.elevation_m=interp2(BATI.x,BATI.y,BATI.h,hazard.lon,hazard.lat);
    
end

% 2) create the storm surge (TS) hazard event set
% -----------------------------------------------

% start from the 'mother' TC hazard event set

if isfield(hazard,'onLand'),hazard=rmfield(hazard,'onLand');end % re-create to be on the safe side
hazard.onLand=hazard.lon.*0+1; % allocate
hazard.onLand(hazard.elevation_m<0)=0; % water points

hazard.peril_ID='TS'; % replace TC with TS
if isfield(hazard,'windfield_comment'),hazard=rmfield(hazard,'windfield_comment');end % remove
if exist('BATI','var')
    hazard.surgefield_comment=sprintf('created based on TC using proxy surge height and bathymetry %s',BATI.sourcefile);
else
    hazard.surgefield_comment=sprintf('created based on TC using proxy surge height and elevation_m');
end
hazard.comment=sprintf('TS hazard event set, generated %s',datestr(now));

% map windspeed onto surge height (the CORE_CONVERSION, see at botton of file, too)
% ===============================
arr_nonzero=find(hazard.intensity); % to avoid de-sparsify all elements
hazard.intensity(arr_nonzero)=0.1023*(max(hazard.intensity(arr_nonzero)-26.8224,0))+1.8288; % m/s converted to m surge height

hazard.elevation_m=max(hazard.elevation_m,0); % only points above sea level

% subtract elevation above sea level from surge height
t0       = clock;
n_events=size(hazard.intensity,1);
n_centroids=size(hazard.intensity,2);

% as the innermost loop is vectorized, it shall be the one repeated most:

if n_events<n_centroids % loop over events, since less events than centroids
    
    msgstr   = sprintf('processing %i events',n_events);
    if climada_global.waitbar
        fprintf('%s (updating waitbar with estimation of time remaining every 100th event)\n',msgstr);
        h        = waitbar(0,msgstr);
        set(h,'Name','Hazard TS: tropical cyclones surge');
    else
        fprintf('%s (waitbar suppressed)\n',msgstr);
        format_str='%s';
    end
    mod_step = 10; % first time estimate after 10 tracks, then every 100
    
    for event_i=1:n_events
        arr_i=find(hazard.intensity(event_i,:)); % to avoid de-sparsify all elements
        hazard.intensity(event_i,arr_i)=max(hazard.intensity(event_i,arr_i)-double(hazard.elevation_m(arr_i)),0); % 20160516 double(.) 
        if mod(event_i,mod_step)==0
            mod_step = 100;
            t_elapsed = etime(clock,t0)/event_i;
            n_remaining = n_events-event_i;
            t_projected_sec = t_elapsed*n_remaining;
            if t_projected_sec<60
                msgstr = sprintf('est. %3.0f sec left (%i/%i events)',t_projected_sec, event_i, n_events);
            else
                msgstr = sprintf('est. %3.1f min left (%i/%i events)',t_projected_sec/60, event_i, n_events);
            end
            
            if climada_global.waitbar
                waitbar(event_i/n_events,h,msgstr); % update waitbar
            else
                fprintf(format_str,msgstr);
                format_str=[repmat('\b',1,length(msgstr)) '%s'];
            end
            
        end
        
    end % event_i
    
else % loop over centroids, since less centroids than events
    
    msgstr   = sprintf('processing %i centroids',n_centroids);
    if climada_global.waitbar
        fprintf('%s (updating waitbar with estimation of time remaining every 100th centroid)\n',msgstr);
        h        = waitbar(0,msgstr);
        set(h,'Name','Hazard TS: tropical cyclones surge');
    else
        fprintf('%s (waitbar suppressed)\n',msgstr);
        format_str='%s';
    end
    mod_step = 10; % first time estimate after 10 tracks, then every 100
    
    for centroid_i=1:n_centroids
        arr_i=find(hazard.intensity(:,centroid_i)); % to avoid de-sparsify all elements
        hazard.intensity(arr_i,centroid_i)=max(hazard.intensity(arr_i,centroid_i)-hazard.elevation_m(centroid_i),0);
        
        if mod(centroid_i,mod_step)==0
            mod_step = 100;
            t_elapsed = etime(clock,t0)/centroid_i;
            n_remaining = n_centroids-centroid_i;
            t_projected_sec = t_elapsed*n_remaining;
            if t_projected_sec<60
                msgstr = sprintf('est. %3.0f sec left (%i/%i centroids)',t_projected_sec, centroid_i, n_centroids);
            else
                msgstr = sprintf('est. %3.1f min left (%i/%i centroids)',t_projected_sec/60, centroid_i, n_centroids);
            end
            
            if climada_global.waitbar
                waitbar(centroid_i/n_centroids,h,msgstr); % update waitbar
            else
                fprintf(format_str,msgstr);
                format_str=[repmat('\b',1,length(msgstr)) '%s'];
            end
        end
        
    end % event_i
    
end
if climada_global.waitbar
    close(h) % dispose waitbar
else
    fprintf(format_str,''); % move carriage to begin of line
end

t_elapsed = etime(clock,t0);
msgstr    = sprintf('generating %i surge fields took %3.2f min (%3.2f sec/event)',n_events,t_elapsed/60,t_elapsed/n_events);
fprintf('%s\n',msgstr);
hazard.creation_comment = msgstr;

if isfield(hazard,'filename'),hazard.filename_source=hazard.filename;end
hazard.filename=hazard_set_file;
hazard.date=datestr(now);
hazard.matrix_density=nnz(hazard.intensity)/numel(hazard.intensity);
hazard.units='m'; % store the SI unit of the hazard intensity
if ~isfield(hazard,'orig_event_count') % fix a minor issue with some hazard sets
    if isfield(hazard,'orig_event_flag')
        fprintf('field hazard.orig_event_count inferred from hazard.orig_event_flag\n')
        hazard.orig_event_count=sum(hazard.orig_event_flag);
    else
        fprintf('WARNING: no field hazard.orig_event_flag\n')
    end
end

if isempty(strfind(hazard_set_file,'NO_SAVE'))
    fprintf('saving TS surge hazard set as %s\n',hazard_set_file);
    save(hazard_set_file,'hazard',climada_global.save_file_version);
end

%%fprintf('TS: max(max(hazard.intensity))=%f\n',full(max(max(hazard.intensity)))); % a kind of easy check

if check_plots,figure;climada_hazard_plot(hazard,0);end % show max surge over ALL events

end % climada_ts_hazard_set

% % =====================================
% % wind speed to surge height conversion (CORE_CONVERSION)
% % =====================================
% 
% % uncomment one level, then copy/paste below into MATLAB command window to run
% 
% mph2ms=0.44704;
% f2m=0.3048;
% 
% % the points read from the SLOSH graph
% v0=60*mph2ms;
% v1=140*mph2ms;
% s0=6*f2m;
% s1=18*f2m;
% 
% % the parameters for the linear function
% a=(s1-s0)/(v1-v0)
% s0
% v0
% 
% figure('Name','windspeed to surge height conversion','Color',[1 1 1])
% hold on
% v=20:100;       plot(v,          a*(v-v0)          +s0,'-r','LineWidth' ,3);
% vmph=60:20:140; plot(vmph*mph2ms,a*(vmph*mph2ms-v0)+s0,'.b','MarkerSize',10);
% legend('conversion','SLOSH points')
% xlabel('wind speed [m/s]')
% ylabel('surge height [m]')
