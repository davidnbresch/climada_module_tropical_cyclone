function hazard = climada_hazard_climate_screw(hazard, hazard_set_file, reference_year, screw)
% NAME:
% climada_hazard_climate_screw
% PURPOSE:
%   Compute the future hazard set for a given climate change scenario from
%   a hazard set defining the risk today.
% CALLING SEQUENCE:
%   hazard = climada_hazard_climate_screw(hazard, hazard_set_file, reference_year, screw)
% EXAMPLE:
%   hazard_cc = climada_hazard_climate_screw(hazard, '', 2040, screw)
%   hazard_cc = climada_hazard_climate_screw(hazard, 'NO_SAVE')
%   hazard_cc = climada_hazard_climate_screw
% INPUTS:
%   hazard:     hazard set today (can be any peril which may be subject to
%               climate change)
% OPTIONAL INPUT PARAMETERS:
%   hazard_set_file:    file name defining where the new hazard should be
%                       saved. If set to 'NO_SAVE', the hazard set will not 
%                       be saved. Prompted for if left empty.
%   reference_year:     time horizon for which one wishes to compute the
%                       climate change hazard set. If left empty, the 
%                       reference year in climada_global will be used.
%   screw:      defines the climate change scenario. A 1xN structure with
%               fields:
%                   .hazard_fld     defines the hazard field to be changed
%                   .change         extent of the change at time horizon
%                   .year           time horizon
%                   .hazard_crit    hazard field to which criteria apply
%                   .criteria       criteria for events/locations to change
%                   .bsxfun_op      operation of change (e.g. @times,@plus) (function handle)
%               specifying N transformations to the original hazard set.
% OUTPUTS:
%   hazard:     new hazard for a future reference year, given a climate
%               change scenario specified by screw.
% MODIFICATION HISTORY:
% Gilles Stassen, gillesstassen@hotmail.com, 20150421 based on original
%           function climada_hazard_clim_scen_advanced by David N. Bresch & Lea Mueller
% Lea Mueller, muellele@gmail.com, 20151021, do not change hazard if it corresponds already to the required year
% Samuel Eberenz, samweli@posteo.de, 20171109, debugged for func2str returning string without leading '@'.
%-

global climada_global
if ~climada_init_vars, return; end

% poor man's version to check arguments
if ~exist('hazard'          , 'var'), hazard           = []; end
if ~exist('hazard_set_file' , 'var'), hazard_set_file  = []; end
if ~exist('reference_year'  , 'var'), reference_year   = []; end
if ~exist('screw'           , 'var'), screw            = []; end

% prompt for hazard if not given
if isempty(hazard) % local GUI
    hazard               = [climada_global.data_dir filesep 'hazards' filesep '*.mat'];
    hazard_default       = [climada_global.data_dir filesep 'hazards' filesep 'choose a hazard.mat'];
    [filename, pathname] = uigetfile(hazard, 'Open existing hazard event set:',hazard_default);
    if isequal(filename,0) || isequal(pathname,0)
        return; % cancel
    else
        hazard = fullfile(pathname,filename);
    end
end
% load the hazard, if a filename has been passed
if ~isstruct(hazard)
    hazard_file = hazard;
    hazard      = [];
    vars = whos('-file', hazard_file);
    load(hazard_file);
    if ~strcmp(vars.name,'hazard')
        hazard = eval(vars.name);
        clear (vars.name)
    end
end

if ~ismember(hazard.peril_ID,{'TC' 'RF' 'TR' 'MA' 'LS'})
    cprintf([1 0 0], sprintf('ERROR: climate change does not apply to %s peril \n',hazard.peril_ID));
elseif strcmp(hazard.peril_ID,'FL')
    fprintf('\t consider running climada_hazard_climate_screw on an RF or TR hazard \n')
    fprintf('set, and computing the screwed FL hazard set from this\n')
    fprintf('aborting\n')
end

%check for statistics and delete
stat_flds   = {'arr_sort' 'arr_ori_sort' 'R' 'R_ori' 'intensity_fit_ori' 'intensity_fit' 'R_fit'};
stat_ndx    = isfield(hazard, stat_flds);
if any(stat_ndx),   hazard = rmfield(hazard, stat_flds(stat_ndx));      end

% use default impact on TC if no screw input given
if ~isstruct(screw) && strcmp(hazard.peril_ID,'TC')
    %enhance frequency of category 4 and 5 storms
    screw(1).hazard_fld         = 'frequency';
    screw(1).change             = 1.8;
    screw(1).year               = 2100;
    screw(1).hazard_crit        = 'category';
    screw(1).criteria           = [4 5];
    screw(1).bsxfun_op          = @times;

    %decrease overall frequency of all storms (0 up to category 5)
    screw(2).hazard_fld         = 'frequency';
    screw(2).change             = 0.72;
    screw(2).year               = 2100;
    screw(2).hazard_crit        = 'category';
    screw(2).criteria           = [0 1 2 3];
    screw(2).bsxfun_op          = @times;
    
    %increase global mean maximum wind speed (intensity)
    screw(3).hazard_fld         = 'intensity';
    screw(3).change             = 1.11;
    screw(3).year               = 2100;
    screw(3).hazard_crit        = 'category';
    screw(3).criteria           = [0 1 2 3 4 5];
    screw(3).bsxfun_op          = @times;
    
    cprintf([0 0 1], 'NOTE: using default (IPCC) values for climate change impact on TC\n')
elseif ~isstruct(screw) && ~strcmp(hazard.peril_ID,'TC')
    cprintf([1 0 0], 'ERROR: climate change screw required as input for non-TC hazard sets\n')
    return;
end

if isempty(reference_year)
    reference_year = climada_global.future_reference_year;
end
fprintf('Reference year for hazard_cc: %d \n',reference_year)

% copy hazard
hazard_cc = hazard;

% implement changes of climate scenario
for i = 1:length(screw)
    % identify relevant events/centroids to change
    crit_ndx = ismember(hazard.(screw(i).hazard_crit),screw(i).criteria);
    
    % linearly interpolate from screw.year to desired reference year
    if (reference_year - hazard.reference_year)~= 0
        time_frac = (reference_year-hazard.reference_year)/(screw(i).year-hazard.reference_year);
        switch func2str(screw(i).bsxfun_op)
            case 'times'
                change = 1+ (screw(i).change-1) * time_frac;
            case '@times'
                change = 1+ (screw(i).change-1) * time_frac;
            case 'plus'
                change = screw(i).change * time_frac;
            case '@plus'
                change = screw(i).change * time_frac;
            otherwise 
                change = screw(i).change;
                warning('screw(%i).bsxfun_op not recognized.',i)
        end % crew(i).bsxfun_op
        fprintf('Apply change %2.4f \n', change)
        
        % determine to which index crit_ndx corresponds (particularly relevant
        % for intensity field), and change by desired amount
        [fld_size_x, fld_size_y] = size(hazard_cc.(screw(i).hazard_fld));
        if length(crit_ndx) == fld_size_x
            hazard_cc.(screw(i).hazard_fld)(crit_ndx,:) = bsxfun(screw(i).bsxfun_op, ...
                hazard.(screw(i).hazard_fld)(crit_ndx,:),change);
        elseif length(crit_ndx) == fld_size_y
            hazard_cc.(screw(i).hazard_fld)(:,crit_ndx) = bsxfun(screw(i).bsxfun_op, ...
                hazard.(screw(i).hazard_fld)(:,crit_ndx),change);
        else
            % something went wrong with indexing...
            cprintf([1 0 0],'ERROR: something went wrong... \n')
            return
        end % length(crit_ndx) == fld_size_x
    else
        fprintf('Hazard corresponds already to request reference year (%d)\n', reference_year)
    end

end % i = 1:length(screw)

% change reference_year and filename
hazard_cc.reference_year = reference_year;
hazard_cc.filename       = [hazard.filename '_cc_' int2str(reference_year)];
hazard_cc.comment        = [hazard.peril_ID ' climate change scenario ' int2str(reference_year)];
hazard                   = hazard_cc;

% prompt for where to save hazard_clim_file if not given
if isempty(hazard_set_file) % local GUI
    hazard_set_file = [climada_global.data_dir filesep 'hazards' filesep '*.mat'];
    if ~exist('filename','var'); filename = '_clim'; else filename = [strtok(filename,'.') '_clim'];end
    hazard_clim_default  = [climada_global.data_dir filesep 'hazards' filesep 'Save climate change hazard in ' filename '.mat'];
    [filename, pathname] = uiputfile(hazard_set_file, 'Save climate change scenario hazard event set as:',hazard_clim_default);
    if isequal(filename,0) || isequal(pathname,0)
        return; % cancel
    else
        hazard_set_file = fullfile(pathname,filename);
    end
    save(hazard_set_file,'hazard')
    cprintf([113 198 113]/255, 'climate change scenario hazard set saved in %s\n', [hazard_set_file]);
elseif ~strcmp(hazard_set_file,'NO_SAVE')
    save(hazard_set_file,'hazard')
    cprintf([113 198 113]/255, 'climate change scenario hazard set saved in %s\n', [hazard_set_file]);
end

end


