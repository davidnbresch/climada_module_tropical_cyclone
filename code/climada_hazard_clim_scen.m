
function hazard = climada_hazard_clim_scen(hazard, tc_track, hazard_save_name, reference_year, screw)

% NAME:
%   climada_hazard_clim_scen
% PURPOSE:
%   starting from a given hazard event set (hazard), construct the
%   hazard event set (hazard_clim_file) for a climate change scenario, e.g.
%   from IPCC SREX
% CALLING SEQUENCE:
%   hazard = climada_tc_hazard_clim_scen(hazard, tc_track, reference_year, screw)
% EXAMPLE:
%   hazard = climada_tc_hazard_clim_scen
% INPUTS:
%   hazard  :         either a hazard set (struct) or a hazard set file (.mat with a struct)
%                     > promted for if not given
%   tc_track:         either a tc track set (struct) or a tc track set file (.mat with a struct)
%                     > promted for if not given (important for storm
%                     category)
% OPTIONAL INPUT PARAMETERS:
%   reference_year:   the reference year for the give climate change scenario, 
%                     e.g. 2017 (+5 years), and then frequency and or intensity 
%                     changes are linearly interpolated from the projected time 
%                     horizon (e.g. 2100 from IPCC SREX) to the requested
%                     reference year.
%   screw:            structure of one or multiple frequency or intensity 
%                     changes given a certain projected time horizon, default
%                     from IPCC SREX taken if not given
%                     example:  
%                        screw.variable_to_change = 'frequency';
%                        screw.frequency          = 0.8;
%                        screw.time_horizon       = 2100;
%                        screw.cat                = [4 5];
%   hazard_save_name: the filename of the new climate scenario hazard event set
%                     > promted for if not given
% OUTPUTS:
%   hazard  :        the hazard event set for the climate scenario, also
%                     stored to hazard_save_name
% MODIFICATION HISTORY:
% David N. Bresch, david.bresch@gmail.com, 20090920
% Lea Mueller, 20120816
%-

global climada_global
if ~climada_init_vars, return; end

% poor man's version to check arguments
if ~exist('hazard'          , 'var'), hazard           = []; end
if ~exist('tc_track'        , 'var'), tc_track         = []; end
if ~exist('hazard_save_name', 'var'), hazard_save_name = []; end
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
        hazard_R_file = fullfile(pathname,filename);
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

%check for statistics and delete
if isfield(hazard,'arr_sort'         ); hazard = rmfield(hazard, 'arr_sort'         );end
if isfield(hazard,'arr_ori_sort'     ); hazard = rmfield(hazard, 'arr_ori_sort'     );end
if isfield(hazard,'R'                ); hazard = rmfield(hazard, 'R'                );end
if isfield(hazard,'R_ori'            ); hazard = rmfield(hazard, 'R_ori'            );end
if isfield(hazard,'intensity_fit_ori'); hazard = rmfield(hazard, 'intensity_fit_ori');end
if isfield(hazard,'intensity_fit'    ); hazard = rmfield(hazard, 'intensity_fit'    );end
if isfield(hazard,'R_fit'            ); hazard = rmfield(hazard, 'R_fit'            );end


% prompt for tc_track if not given
if isempty(tc_track) % local GUI
    tc_track               = [climada_global.data_dir filesep 'tc_tracks' filesep '*.mat'];
    tc_track_default       = [climada_global.data_dir filesep 'tc_tracks' filesep 'choose a tc_track.mat'];
    [filename, pathname]   = uigetfile(tc_track, 'Open existing tc track set:',tc_track_default);
    if isequal(filename,0) || isequal(pathname,0)
        return; % cancel
    else
        tc_track = fullfile(pathname,filename);
    end
end
% load the tc_track, if a filename has been passed
if ~isstruct(tc_track)
    tc_track_file = tc_track;
    tc_track      = [];
    vars = whos('-file', tc_track_file);
    load(tc_track_file);
    if ~strcmp(vars.name,'tc_track')
        tc_track = eval(vars.name);
        clear (vars.name)
    end
end
% add tc track category
if ~isfield(tc_track, 'category')
    tc_track = climada_tc_stormcategory(tc_track);
    fprintf('field tc_track.category added \n')
end

if isempty(screw)
    %enhance frequency of category 4 and 5 storms
    screw.variable_to_change = 'frequency';
    screw.frequency          = 0.8;
    screw.time_horizon       = 2100;
    screw.cat                = [4 5];
    
    %decrease overall frequency of all storms (0 up to category 5)
    screw(2).variable_to_change = 'frequency';
    screw(2).frequency          = -0.28;
    screw(2).time_horizon       = 2100;
    screw(2).cat                = [0 1 2 3];
    
    %increase global mean maximum wind speed (intensity)
    screw(3).variable_to_change = 'intensity';
    screw(3).frequency          = 0.11;
    screw(3).time_horizon       = 2100;
    screw(3).cat                = [0 1 2 3 4 5];
    
    fprintf('***\nDefault values for climate change:\n')
    %fprintf('screw.frequency   : %10.2f%%\n',screw(1).frequency*100)
    %fprintf('screw.time_horizon: %d      \n',screw(1).time_horizon)
    %fprintf('screw.cat         : %s      \n',int2str(screw(1).cat))
end


if isempty(reference_year)
    reference_year     = 2017;
end
fprintf('Reference year for hazard_cc: %d \n',reference_year)


%% copy hazard
hazard_cc = hazard;


 

%% implement changes of climate scenario
for screw_i = 1:length(screw)
    frequency_screw_reference_year ...
                   = (reference_year              - hazard.reference_year)...
                    /(screw(screw_i).time_horizon - hazard.reference_year)...
                    * screw(screw_i).frequency;
    
    % put all categories of the tc tracks in one vector
    tc_all_cat       = [tc_track.category];
    tc_to_change     = ismember(tc_all_cat,screw(screw_i).cat);

    if strcmp(screw(screw_i).variable_to_change,'frequency')
        % change the frequency
        hazard_cc.frequency(tc_to_change) = hazard_cc.frequency(tc_to_change)*(1+frequency_screw_reference_year);
        
    elseif strcmp(screw(screw_i).variable_to_change,'intensity')
        % change the intensity
        hazard_cc.arr(tc_to_change,:) = hazard_cc.arr(tc_to_change,:)*(1+frequency_screw_reference_year);
    end
    
    if frequency_screw_reference_year>0
        change_action = 'increased';
    else
        change_action = 'decreased';
    end
    fprintf('***\n%s %s by %10.2f%% for category %s for reference year %d\n',...
                             screw(screw_i).variable_to_change,...
                             change_action, ...
                             frequency_screw_reference_year*100, ...
                             int2str(screw(screw_i).cat),...
                             reference_year)           
end


                     
%% change reference_year and filename
hazard_cc.reference_year = reference_year;
hazard_cc.filename       = [hazard.filename '_cc_' int2str(reference_year)];
hazard_cc.comment        = ['TCNA climate change scenario ' int2str(reference_year)];
hazard                   = hazard_cc; 

% prompt for where to save hazard_clim_file if not given
if isempty(hazard_save_name) % local GUI
    hazard_save_name = [climada_global.data_dir filesep 'hazards' filesep '*.mat'];
    if ~exist('filename','var'); filename = '_clim'; else filename = [strtok(filename,'.') '_clim'];end
    hazard_clim_default  = [climada_global.data_dir filesep 'hazards' filesep 'Save climate change hazard in ' filename '.mat'];
    [filename, pathname] = uiputfile(hazard_save_name, 'Save climate change scenario hazard event set as:',hazard_clim_default);
    if isequal(filename,0) || isequal(pathname,0)
        return; % cancel
    else
        hazard_save_name = fullfile(pathname,filename);
    end
else
    hazard_save_name = [climada_global.data_dir filesep 'hazards' filesep hazard_save_name];
end
save(hazard_save_name,'hazard')
fprintf('\n***Climate change scenario *** \n');
cprintf([113 198 113]/255, 'saved in %s\n', [hazard_save_name]);




