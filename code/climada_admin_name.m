function [admin_name, admin_shapes, adm_lvl, shape_ndx, country_name, location] ...
    = climada_admin_name(country_name,admin_name,adm_lvl,climada_global_border_file_check)
% climada_admin_name
% MODULE:
%   climada core
% NAME:
%   climada_admin_name
% PURPOSE:
%   check for valid admin name on different levels (0: country, 1: states, 2: districts) 
%   and download admin shape files, return all the arguments incl. shape files
%   Same concept as climada_country_name, but for admin regions. The map
%   border file in climada_global as used in climada_plot_world_borders
%   may be replaced with the shape files for the desired admin level, but
%   this may cause problems later
% CALLING SEQUENCE:
%   [admin_name, admin_shapes, adm_lvl, shape_ndx, country_name, location] ...
%     = climada_admin_name(country_name,admin_name,adm_lvl,climada_global_border_file_check)
% EXAMPLE:
%   [admin_name, admin_shapes, adm_lvl, shape_ndx, country_name, location] ...
%     = climada_admin_name('Netherlands','Utrecht',2,1)
% INPUTS: 
% OPTIONAL INPUT PARAMETERS:
%   country_name:   country of interest,prompted from if not given
%   admin_name:     set to 'SINGLE' (default) for single admin region
%                   selection at the specified admin level
%                   set to 'MULTIPLE' for multi-region select within one
%                   admin level, e.g. admin level set to 2 => select
%                   country, select admin region at level 1, then choose
%                   multiple regions at level 2
%                   set to 'ALL' to return all admin regions at desired
%                   admin level for entire country
%   adm_lvl:        define admin level of interest (max admin level varies
%                   per country)
%   climada_global_border_file_check:   whether to set admin shapes as
%                                       default climada global map border file.
% OUTPUTS:
%   admin_name:     name of the selected admin region
%   admin_shapes:   shapes struct array at admin level
%   adm_lvl:        administrative region (0 = country)
%   shape_ndx:      shape index of selected admin region in admin shapes
%                   structure array
%   country_name:   country name
%   location:       struct with .lon .lat for centre of selected admin region
% MODIFICATION HISTORY:
% Gilles Stassen, gillesstassen@hotmail.com, 20150508 init
% Gilles Stassen, gillesstassen@hotmail.com, 20150615, documentation, hopefully all bugs are gone (they never are...)
% Lea Mueller, muellele@gmail.com, 20150615, climada_global.data dir instead of module_dir
%-
admin_shapes = []; shape_ndx = []; location = []; % init output

global climada_global
if ~climada_init_vars,return;end % init/import global variables

% need correct border file for climada country name
c_g_switch = 0;
if isfield(climada_global,'climada_global_ori')
    c_g_switch = 1;
    tmp_climada_global = climada_global;
    climada_global = climada_global.climada_global_ori;
end

% Check input variables
if ~exist('country_name','var') || isempty(country_name)
    [country_name, ISO3, country_shape_ndx] = climada_country_name('SINGLE');
    if isempty(country_name), admin_name = []; return; end
    country_name    = char(country_name);
    ISO3            = char(ISO3);
else
    [country_name, ISO3, country_shape_ndx] = climada_country_name(country_name);
    if isempty(country_name)
        cprintf([1 0 0],'ERROR: invalid country name \n');
        admin_name = [];
        return; 
    end
end

if c_g_switch
    climada_global = tmp_climada_global; clear tmp_climada_global
end

if ~exist('admin_name',             'var'),	admin_name = 'SINGLE';      end
if ~exist('adm_lvl',                'var'),	adm_lvl='';               	end
if ~exist('climada_global_border_file_check','var'), climada_global_border_file_check = 0;  end

if climada_global_border_file_check
    cprintf([0 0 1], ['NOTE: \tchanging the map border file stored in climada_global '...
        'may cause issues \n\t\twith some functions, including climada_country_name\n'])
    if ~isfield(climada_global,'climada_global_ori')
        climada_global.climada_global_ori = climada_global;
    end
end

% the shape file with higher resolution (than default climada) for country
country_shapefile = [climada_global.data_dir filesep 'entities' filesep ISO3 '_adm' filesep ISO3 '_adm0.shp'];

% download border shapefiles
admin_dir = fileparts(country_shapefile);
if ~exist(country_shapefile,'file')
    admin_URL = ['http://biogeo.ucdavis.edu/data/diva/adm/' ISO3 '_adm.zip'];
    fprintf('downloading admin region shapefiles from http://www.diva-gis.org/gdata... ');
    unzip(admin_URL,admin_dir);
    fprintf('done \n')
end

files = dir(admin_dir);
adm_max = -1; % init (adm lvls start at 0)
for f_i = 1: length(files)
    [~,~,fE] = fileparts([admin_dir filesep files(f_i).name]);
    % still need these for field names and attributes!!
    %     if ~(strcmp(fE,'.mat') || strcmp(fE,'.shp') || strcmp(fE,'.'))
    %         delete([admin_dir filesep files(f_i).name]);
    %     end
    if strcmp(fE,'.shp');   adm_max = adm_max +1;     end
end

% Specify the administrative level of interest (from function input)
% The data was originally downloaded from http://www.diva-gis.org/gdata
if isempty(adm_lvl)
    adm_lvl = input(sprintf('Choose the admin level of interest (0-%i): ',adm_max));
end
if adm_lvl > adm_max || adm_lvl < 0
    fprintf('ERROR: please choose an admin level from 0 to %i, aborted',adm_max);
    return;
end
% redefine filename according to chosen admin level
admin_regions_shapefile = [climada_global.data_dir filesep 'entities' filesep ISO3 '_adm' filesep ISO3 '_adm' num2str(adm_lvl) '.shp'];
clear adm_count files f_i fE

% If the file specified above exists, use the high resolution borders
% stored there, else revert to default climada low resolution borders.
if exist(admin_regions_shapefile,'file')
    % admin_shapes = shaperead(admin_regions_shapefile);
    admin_shapes = climada_shaperead(admin_regions_shapefile,1);
    
    fld_ID      = sprintf('ID_%i',      adm_lvl);
    fld_NAME    = sprintf('NAME_%i',    adm_lvl);
    
    if (isempty(admin_name) || strcmp(admin_name,'SINGLE') || strcmp(admin_name,'MULTIPLE'))...
            && adm_lvl ~=0
%--------------------------------        
        shape_ndx_s = 1:length(admin_shapes);
        chosen_admin={}; % init
        for adm_i=1:adm_lvl
            % select from map
            admin_regions_shapefile = [climada_global.data_dir filesep 'entities' filesep ISO3 '_adm' filesep ISO3 '_adm' num2str(adm_i) '.shp'];
            admin_shapes = climada_shaperead(admin_regions_shapefile,1);
            shape_ndx_s = 1:length(admin_shapes);
            if ~isempty(chosen_admin)
                for shape_i = 1: length(admin_shapes)
                    if ~strcmp(admin_shapes(shape_i).(fld_NAME_tmp),chosen_admin(adm_i-1))
                        shape_ndx_s(shape_i) = NaN;
                    end
                end
                shape_ndx_s = shape_ndx_s(~isnan(shape_ndx_s));
                admin_shapes = admin_shapes(shape_ndx_s);
            end
            fld_NAME_tmp    = sprintf('NAME_%i',    adm_i);
            
            msg_str = sprintf('select admin level %i region of interest',adm_i);
            figure('name',msg_str,'color','w','outerposition',[200 148 1000 720]); hold on
%             axis([shapes(country_shape_ndx).BoundingBox(:,1)' shapes(country_shape_ndx).BoundingBox(:,2)'])
            for shape_i = 1: length(admin_shapes)
                % plot borders
                plot(admin_shapes(shape_i).X,admin_shapes(shape_i).Y,'color',[81 81 81]./255,'linewidth',2)
                % label admin regions
                NAME_lon = mean(admin_shapes(shape_i).BoundingBox(:,1));
                NAME_lat = mean(admin_shapes(shape_i).BoundingBox(:,2));
                t = text(NAME_lon,NAME_lat,admin_shapes(shape_i).(fld_NAME_tmp));
                set(t,'HorizontalAlignment','center','color','b');
            end
            axis tight
            axis equal
            axis off
            if strcmp(admin_name,'MULTIPLE') && adm_i == adm_lvl
                fprintf('click on regions of interest (press enter when done): ')
                [click_lon,click_lat] = ginput;
            else
                fprintf('click on region of interest: ')
                [click_lon,click_lat] = ginput(1);
            end
            close
            for click_i = 1:length(click_lon)
                for shape_i = 1: length(admin_shapes)
                    if inpolygon(click_lon(click_i),click_lat(click_i),...
                            admin_shapes(shape_i).X,admin_shapes(shape_i).Y)
                        chosen_admin{adm_i+click_i-1} = admin_shapes(shape_i).(fld_NAME_tmp);
                        break;
                    end
                end
                fprintf(['\t' chosen_admin{adm_i+click_i-1} ]);
            end
            fprintf('\n')
%--------------------------------
% paste dialogue box routine here, and remove code above - not sure if it
% still works though....
%--------------------------------
        end
        admin_name = chosen_admin(adm_lvl:end);
        
        % index for label
        for click_i = 1:length(click_lon)
            for shape_i = 1 : length(admin_shapes)
                if strcmp(admin_shapes(shape_i).(fld_NAME), admin_name(click_i))  %strfind(eval(strcat('admin_shapes(i).NAME_',num2str(adm_lvl))), admin_name)
                    ID = admin_shapes(shape_i).(fld_ID);  %ID = eval(strcat('admin_shapes(i).ID_',num2str(adm_lvl)));
                    shape_ndx(click_i) = shape_i;
                    break;
                end
                if shape_i == numel(admin_shapes), fprintf('ERROR: %s not found',admin_name); return; end
            end
        end
    elseif (isempty(admin_name) && adm_lvl == 0) || strcmp(admin_name,'ALL')
        admin_name = country_name;
    else
        for shape_i = 1 : length(admin_shapes)
            if strcmp(admin_shapes(shape_i).(fld_NAME), admin_name)  %strfind(eval(strcat('admin_shapes(i).NAME_',num2str(adm_lvl))), admin_name)
                ID = admin_shapes(shape_i).(fld_ID);  %ID = eval(strcat('admin_shapes(i).ID_',num2str(adm_lvl)));
                shape_ndx = shape_i;
                break;
            end
            if shape_i == numel(admin_shapes), fprintf('ERROR: %s not found',admin_name); return; end
        end
    end
    
    %     clear shape_ndx_s fld_NAME_tmp
    
    % make sure there is a .mat version of the border shapes file
    [fP,fN] = fileparts(country_shapefile);
    country_shapefile_mat=[fP filesep fN '.mat'];
    if ~exist(country_shapefile_mat,'file')
        % first time, read the shape file and store as .mat
        country_shapes = climada_shaperead(country_shapefile);
    else
        load(country_shapefile_mat)
        country_shapes = shapes;
    end
    
    % annotation (just to label admin region on plots):
    if adm_lvl > 0 && exist('shape_ndx','var') && nargout == 4
        if ~iscell(admin_name),     admin_name = {admin_name};  end
        for loc_i = 1:length(admin_name)
            location(loc_i).name   = ['  ' admin_name{loc_i}]; % first two spaces for nicer labeling
            location(loc_i).lon    = median(admin_shapes(shape_ndx(loc_i)).X);
            location(loc_i).lat    = median(admin_shapes(shape_ndx(loc_i)).Y);
        end
    end
    
    % for plotting routines, set admin shapes as default border file
    [fP,fN] = fileparts(admin_regions_shapefile);
    
    if climada_global_border_file_check
        admin_regions_matfile=[fP filesep fN '.mat'];
        climada_global.map_border_file = admin_regions_matfile; % set map border file for nicer plotting 
        try
            figure; climada_plot_world_borders; close
        catch
            cprintf([1 0.5 0],spritnf('WARNING: unable to set %s.mat as border file\n',fN))
        end
    end
    
%     admin_shapes = admin_shapes(shape_ndx);
    
    clear NAME_lat NAME_lon click_lat click_lon shape_ndx_s fld_NAME_tmp bounding_box shape_i

else
    fprintf('WARNING: country shapefiles not found, using standard resolution world map instead\n');
    fprintf('\t download shapefiles from: www.diva-gis.org/gdata \n');
end


%--------------------------------
% the following lines of code allow the user to select a region of interest
% from a list dialoge box (as opposed to clicking on the map) - paste above
%
%             % construct entries for listdlg
%             admin_names={}; % init
%             fld_NAME_tmp    = sprintf('NAME_%i',    adm_i);
%             for shape_i=1:length(admin_shapes)
%                 if ismember(shape_i,shape_ndx_s)
%                     admin_names{shape_i}=admin_shapes(shape_i).(fld_NAME_tmp);
%                 else
%                     admin_names{shape_i}='';
%                 end
%             end % shape_i
%
%             if any(isempty(admin_names))
%                 admin_names{strcmp(admin_names,'[]')} = '';
%             end
%
%             [liststr,sort_ndx] = unique(admin_names); % sort alphabetically
%             selection = listdlg('PromptString',sprintf('Select admin %i region:',adm_i),...
%                 'ListString',liststr,'SelectionMode','SINGLE');
%             pause(0.1)
%
%             if ~isempty(selection)
%                 admin_name = char(admin_names(sort_ndx(selection)));
%                 shape_ndx  = sort_ndx(selection);
%             else
%                 fprintf('NOTE: no region chosen, aborted\n')
%                 return
%             end
%             shape_ndx_s(~strcmp(admin_names,admin_name)) = NaN;
%--------------------------------